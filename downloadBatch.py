import bz2
import gzip
import multiprocessing
import os
import pandas
import shutil
import subprocess
import sys
import time
import xml.etree.ElementTree as etree

pandas.options.display.width=250
"""
Study
Experiment
Run
(Sample)
"""

def splitFastqStream(stream, f1, f2):
    for i, line in enumerate(stream):
        if i % 8 < 4:
            f1.write(line)
        else:
            f2.write(line)



def runCmd(cmd, verbose=True):
    if verbose:
        print "running cmd: '{}'".format(cmd)

    t0 = time.time()
    p = subprocess.Popen(cmd, shell=True)
    result = p.wait()
    t1 = time.time()

    if result != 0:
        raise Exception("Failed to run command '{}' -- error code {}".format(cmd, result))
    if verbose:
        print "  ...time elapsed: {:.1f}s      [[{}]]".format(t1-t0, cmd)

def findone(node, match):
    found = node.findall(match)
    if len(found) < 1:
        raise Exception("None found: {}: {}".format(node, match))
    if len(found) > 1:
        raise Exception("Multiple found, one expected: {}: {} (n={})".format(node, match, len(found)))
    return found[0]


def loadMetadata(accession):
    metadata = subprocess.check_output("esearch -db sra -query {} | efetch".format(accession), shell=True)
    return etree.fromstring(metadata)

def buildTable(accession):
    table = []
    metadata = loadMetadata(accession)

    # experiment
    for package in metadata.iter("EXPERIMENT_PACKAGE"):
        experimentName = findone(package, "EXPERIMENT").attrib["accession"]
        experimentLong = findone(package, "EXPERIMENT/TITLE").text

        submissionName = findone(package, "SUBMISSION").attrib["accession"]
        studyName = findone(package, "STUDY").attrib["accession"]
        studyLong = findone(package, "STUDY/DESCRIPTOR/STUDY_TITLE").text

        sampleName = findone(package, "SAMPLE").attrib["accession"]
        sampleLong = findone(package, "SAMPLE/TITLE").text.replace(" ", "_")

        isPaired = findone(package, "EXPERIMENT").find("DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_LAYOUT/PAIRED")!=None

        for run in findone(package, "RUN_SET").iter("RUN"):
            readCount = int(run.attrib["total_spots"])
            table.append([submissionName, studyName, studyLong, experimentName, experimentLong, sampleName, sampleLong, 
                run.attrib["accession"], isPaired, readCount])

    return pandas.DataFrame(table, columns=["submission", "study", "study_long", "expt", "expt_long", "sample", 
        "sample_long", "run", "paired", "read_count"])

def downloadTable(table, where, filter=None):
    # make directories first, single-threaded
    for i, row in table.iterrows():
        if not os.path.exists(row["sample_long"]):
            os.makedirs(row["sample_long"])

    # now run the downloads 
    if where == "japan":
        # don't do it in parallel, as the japanese mirror has a limit of 3 concurrent
        # processes
        for i, row in table.iterrows():
            downloadRunFromJapan(row)
    elif where == "ebi":
        for i, row in table.iterrows():
            downloadRunFromEBI(row)        
    else:
        # knock yourself out -- each process is going to be horrendously slow 
        # because of the fastq-dump step, so might as well parallelize like crazy
        pool = multiprocessing.Pool(processes=8)
        pool.map_async(downloadRun_wrapper, [row for i, row in table.iterrows()]).get(999999)


def downloadRun_wrapper(row):
    import traceback
    try:
        t0 = time.time()
        downloadRunFromUSA(row)
        t1 = time.time()
        print "Time elapsed to download {}: {}s".format(row["run"], t1-t0)
    except Exception, e:
        print "Exception on row \n'{}':\n'{}'".format(row, traceback.format_exc())
        raise



def _checkFile(path, row, suffix):
    if not os.path.exists(path):
        return False
    print "verifying path:", path

    fileOpenCmd = open
    if suffix == ".bz2":
        fileOpenCmd = bz2.BZ2File
    elif suffix == ".gz":
        fileOpenCmd = gzip.open

    # try:
    if suffix == ".bz2":
        try:
            subprocess.check_call("which lbzip2", shell=True, stdout=subprocess.PIPE)
        except:
            print "\n"
            print "you should install lbzip2!"
            print "\n"
            raise

        p = subprocess.Popen("lbzip2 -d -c {} | wc -l".format(path), shell=True, stdout=subprocess.PIPE)
        count = int(p.stdout.read().strip()) / 4
    elif suffix == ".gz":
        cmd = "gzip -c -d {} | wc -l".format(path)
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        count = int(p.stdout.read().strip()) / 4
    else:
        with fileOpenCmd(path) as f:
            count = 0
            for i, line in enumerate(f):
                if len(line) > 1:
                    count += 1
            count /= 4
    if count == row["read_count"]:
        return True

    print "@"*10
    print row["sample_long"], "count:", count, "expected:", row["read_count"]
    print "@"*10

    return False


def checkRun(row, suffix=".gz"):
    result = True

    if row["paired"]:
        for end in [1,2]:
            curPath = "{}/{}_{}.fastq{}".format(row["sample_long"], row["run"], end, suffix)
            result &= _checkFile(curPath, row, suffix)
    else:
        curPath = "{}/{}.fastq{}".format(row["sample_long"], row["run"], suffix)
        result &= _checkFile(curPath, row, suffix)

    return result


def downloadRunFromEBI(row):
    """ They ship fastq.gz files directly, none of this .sra business """

    run = row["run"]


    if checkRun(row):
        print "Run alread downloaded and contains the correct number of reads"
        return

    ends = [""]
    if row["paired"]:
        ends = ["_1", "_2"]

    for end in ends:
        inpath = "{}{}.fastq.gz".format(run, end)
        outpath = "{}/{}{}.fastq.gz".format(row["sample_long"], run, end)

        if not os.path.exists(outpath):
            # see http://www.ebi.ac.uk/ena/browse/read-download
            dir1 = row["run"][:6]
            if len(row["run"]) == 9:
                dir2 = ""
            elif len(row["run"]) == 10:
                dir2 = "/00{}".format(row["run"][-1])
            elif len(row["run"]) == 11:
                dir2 = "/0{}".format(row["run"][-2:])
            elif len(row["run"]) == 12:
                dir2 = "/{}".format(row["run"][-3:])

            url = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/{}{}/{}/{}".format(
                dir1, dir2, row["run"], inpath)

            cmd1 = "aria2c -x 16 {}".format(url)
            print "Downloading run: {} Size:{:,} reads ".format(row["sample_long"], row["read_count"], "(paired)" if row["paired"] else "")
            runCmd(cmd1)

            print "moving {} to {}".format(inpath, outpath)
            shutil.move(inpath, outpath)
        else:
            print "already downloaded:{}".format(run)

    if not checkRun(row, ".bz2"):
        print "*"*100
        print "Failed to download and expand the following run:"
        print row
        print "*"*100
        print ""

    print "** Run Complete: {} **".format(run)


def downloadRunFromJapan(row):
    """ They ship fastq.bz2 files directly, none of this .sra business

    HOWEVER, they have a limit of 3 concurrent connections from a single IP """
    run = row["run"]


    if checkRun(row):
        print "Run alread downloaded and contains the correct number of reads"
        return

    ends = [""]
    if row["paired"]:
        ends = ["_1", "_2"]

    for end in ends:
        inpath = "{}{}.fastq.bz2".format(run, end)
        outpath = "{}/{}{}.fastq.bz2".format(row["sample_long"], run, end)

        if not os.path.exists(outpath):
            url = "ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/{}/{}/{}/{}".format(
                row["submission"][:6], row["submission"], row["expt"], inpath)
            cmd1 = "aria2c -x 3 {}".format(url)
            print "Downloading run: {} Size:{:,} reads ".format(row["sample_long"], row["read_count"], "(paired)" if row["paired"] else "")
            runCmd(cmd1)

            print "moving {} to {}".format(inpath, outpath)
            shutil.move(inpath, outpath)
        else:
            print "already downloaded:{}".format(run)

    if not checkRun(row, ".bz2"):
        print "*"*100
        print "Failed to download and expand the following run:"
        print row
        print "*"*100
        print ""

    print "** Run Complete: {} **".format(run)


def downloadRunFromUSA(row):
    """ For some reason, converting from .sra to .fastq.gz takes *forever*, so you're best off downloading the 
    fastq's directly from Japan if you can """
    run = row["run"]

    print "starting..."
    if checkRun(row, suffix=".gz"):
        print "Run alread downloaded and contains the correct number of reads"
        return

    if not os.path.exists("{}.sra".format(run)):
        url = "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/{}/{}/{}.sra".format(
            run[:6], run, run)

        # Not sure if there's a limit to concurrent connections to NIH
        # but this can go up to 16 as far as aria2c is concerned
        cmd1 = "aria2c -x 8 {}".format(url)
        print "Downloading run: {} Size:{:,} reads ".format(row["sample_long"], row["read_count"], "(paired)" if row["paired"] else "")
        runCmd(cmd1)
    else:
        print "already downloaded:{}".format(run)

    inpaths = []
    outpaths = []

    if row["paired"]:
        for end in [1,2]:
            inpaths.append("{}_{}.fastq".format(run, end)) # fix if using --gzip
            outpaths.append("{}/{}_{}.fastq.gz".format(row["sample_long"], run, end))
    else:
        inpaths.append("{}.fastq".format(run)) # fix if using --gzip
        outpaths.append("{}/{}.fastq.gz".format(row["sample_long"], run))        

    fastqsExist = True
    for inpath in inpaths:
        if not os.path.exists(inpath):
            print "missing:{}".format(inpath)
            fastqsExist = False

    if not fastqsExist:
        if row["paired"]:
            # -I means add read-end ID to name of read in header line
            # --split-files is required to make two files for paired end data
            cmd2 = "fastq-dump --split-files -I {}".format(row["run"])

            # this (using --gzip) is slow as all hell:
            # cmd2 = "fastq-dump --split-files -I --gzip {}".format(row["run"])
        else:
            cmd2 = "fastq-dump {}".format(row["run"])

        runCmd(cmd2)


    for inpath, outpath in zip(inpaths, outpaths):
        print "compressing {} to {}".format(inpath, outpath)
        cmd3 = "pigz -c {} > {}".format(inpath, outpath)
        runCmd(cmd3)

        # if using 'fastq-dump --gzip' above
        # shutil.move(inpath, outpath)

    if not checkRun(row, suffix=".gz"):
        print "*"*100
        print "Failed to download and expand the following run:"
        print row
        print "*"*100
        print ""

    print "** Run Complete: {} **".format(run)

def main():
    accessions = sys.argv[1:]
    where = "usa"
    if "--japan" in accessions:
        print "downloading from japan mirror"
        where = "japan"
        accessions = [x for x in accessions if x!="--japan"]
    if "--ebi" in accessions:
        print "downloading from european mirror"
        where = "ebi"
        accessions = [x for x in accessions if x!="--ebi"]

    table = None
    for accession in accessions:
        curtable = buildTable(accession)
        if table is None:
            table = curtable
        else:
            table = pandas.concat([table, curtable])

    print table

    def fileFormat(**x):
        sample_long = x["sample_long"].replace(" ", "_")
        return "{}/{}".format(sample_long, x["run"])
    downloadTable(table, where)

if __name__ == '__main__':
    main()




# python ~/projects/sra/downloadBatch.py GSM1155957 GSM1155958 GSM1155959 GSM1155960

""" Approximate timings (SRX437468):

Japan               6.5m    (missing 30k lines)
USA (with pigz)     14m
USA (with --gzip)   22m


ftp.sra.ebi.ac.uk/vol1/fastq/SRR891/SRR891268/SRR891268_1.fastq.gz
ftp.sra.ebi.ac.uk/vol1/fastq/SRR891/SRR891268/SRR891268_1.fastq.gz
"""