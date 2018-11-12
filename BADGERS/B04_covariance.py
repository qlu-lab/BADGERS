#!/usr/bin/env python

from timeit import default_timer as timer


import logging
import numpy
import os
import gzip
import ntpath
import badgers.WeightDBUtilities as WeightDBUtilities
import badgers.PrediXcanFormatUtilities as PrediXcanFormatUtilities
import badgers.ThousandGenomesUtilities as ThousandGenomesUtilities
import badgers.Logging as Logging
import badgers.Utilities as Utilities
import badgers.Formats as Formats



def pathLeaf(path):
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)

class ProcessWeightDB(object):
    def __init__(self, args):
        self.weight_db = pathLeaf(args.weight_db)
        self.db_path = args.weight_db
        self.data_folder = args.input_folder
        self.correlation_output = args.correlation_output
        self.covariance_output = args.covariance_output
        if args.covariance_output is None:
            comp = os.path.splitext(self.weight_db)[0]
            name = comp + ".cov.txt.gz"
            path = os.path.join("intermediate", "cov")
            path = os.path.join(path, name)
            self.covariance_output = path

        self.input_format = args.input_format

        self.found_genes_for_covariance = {}
        self.found_genes_for_correlation = {}

        self.min_maf_filter = float(args.min_maf_filter) if args.min_maf_filter else None
        self.max_maf_filter = float(args.max_maf_filter) if args.max_maf_filter else None

        self.max_snps_in_gene = int(args.max_snps_in_gene) if args.max_snps_in_gene else None

    def run(self):
        start = timer()
        if os.path.exists(self.covariance_output):
            logging.info("%s already exists, delete it if you want it figured out again", self.covariance_output)
            return
        if not self.correlation_output and not self.covariance_output:
            logging.info("Provide --correlation_output or --covariance_output or both")
            return
        logging.info("Loading Weights")
        weight_db_logic = WeightDBUtilities.WeightDBEntryLogic(self.db_path)
        logging.info("Building files")
        head, tail = os.path.split(self.covariance_output)
        file_name = tail.split(".")[0]
        data_location = head+"/"+file_name+"_tn.txt.gz"
        entries,wat = self.buildFiles(weight_db_logic,data_location)
        logging.info("Ran successfully")
        end = timer()
        logging.info("Sucessfully processed metaxcan covariance in %s seconds"%(str(end - start)))
        return entries,wat

    def buildFiles(self, weight_db_logic,data_location):
        do_covariances = True
        self.writeFileHeader(self.covariance_output)

        names = Utilities.dosageNamesFromFolder(self.data_folder)
        wat = []
        for x in xrange(503):
            wat.append(0);
            x+=1
        for name in names:
            seven = timer()
            snps, snps_by_rsid = self.getSNPS(name, weight_db_logic)
            eight = timer()
            one = timer()
            if do_covariances:
                entries,wat_sub = self.addToCovarianceFile(weight_db_logic, name, snps, snps_by_rsid,data_location)
                for i in xrange(503) :
                    wat[i] += float(wat_sub[i])
                    i+=1  
            two = timer()
        return entries,wat
    def writeFileHeader(self,path):
        with gzip.open(path, "ab") as file:
            file.write("TRAIT RSID1 RSID2 VALUE\n")

    def getSNPS(self, name, weight_db_logic):
        dosageLoader = None
        if self.input_format == Formats.IMPUTE:
            dosageLoader = ThousandGenomesUtilities.IMPUTEDosageLoader(self.data_folder, name) #outdated code
        elif self.input_format == Formats.PrediXcan:
            dosageName = Utilities.dosageName(name)
            path = os.path.join(self.data_folder, dosageName)
            dosageLoader = PrediXcanFormatUtilities.PrediXcanFormatDosageLoader(path, weight_db_logic)
        else:
            logging.info("Invalid input format: %s", self.input_format)
            return
        snps, snps_by_rsid = dosageLoader.load()
        return snps, snps_by_rsid

    def addToCovarianceFile(self, weight_db_logic, name, snps, snps_by_rsid,data_location):
        logging.info("Adding to covariance for %s-%s", name, self.weight_db)
        ka = timer()
        genes = weight_db_logic.weights_by_gene.keys()
        total_genes = len(genes)
        last_reported_percent = 0
        processed = 0
        kb = timer()
        entries = []
        wat = []
        for x in xrange(503):
            wat.append(0);
            x+=1
        for gene in genes:
            processed += 1
            percent = int(processed*100.0 / total_genes)
            if percent == last_reported_percent+10:
                logging.info("%d percent genes processed", percent)
                last_reported_percent = percent
            three=timer()
            entries_sub,wat_sub = self.buildCovarianceEntries(name, gene, weight_db_logic, snps_by_rsid,data_location)
            if len(entries_sub) == 0:
                continue  
            entries = entries + entries_sub
            for i in xrange(503) :
                wat[i] += float(wat_sub[i])
                i+=1  
            four = timer()
            if len(entries) == 0:
                logging.log(6,"Gene %s has no snps in current file", gene)
                continue
            kc = timer()
            self.addToFile(self.covariance_output, gene, entries)
            kd = timer()

        return entries,wat
    def addToFile(self, path, gene, entries):
        get = path.split("/")[-1]
        real = (get[0:-7])
        with gzip.open(path, "ab") as file:
            for entry in entries:
                line = " ".join([real, entry[0], entry[1], entry[2]])+"\n"
                file.write(line)
    def buildCovarianceEntries(self, name, gene, weight_db_logic, snps_by_rsid,data_location):
        weights_in_gene = weight_db_logic.weights_by_gene[gene]
        rsids_from_genes = weights_in_gene.keys()
        nine = timer()
        related_rsids, related_data = self.buildRelatedData(rsids_from_genes, snps_by_rsid, weights_in_gene)
        weight=[]
        for k in related_rsids:
            weight.append(weights_in_gene[k].weight)
        ten = timer()
        if len(related_rsids) == 0:
            return [],[]
        self.updateFoundCovariance(gene, name)
        array = numpy.array(related_data)
        eleven = timer()
        var_mat = []
        for lines in array:
            varianc = numpy.var(lines);
            var_mat.append(varianc)
        var_mat_true = numpy.array(var_mat)

        fw = gzip.open(data_location,'a+')
        for kx in xrange(len(related_rsids)):
            out = [str(i) for i in related_data[kx]]
            data = " ".join(out)
            total = str(related_rsids[kx])+" "+data
            fw.write(total+"\n")
        twelfe = timer()
        ti = numpy.dot(array.transpose(),weight)
        wat = []
        for abc in xrange(len(ti)):        
            wat.append(str(ti[abc]))
            abc+=1
        #translate into sql entries
        five = timer()
        entries = self.buildMatrixOutputEntries(var_mat_true, rsids_from_genes, related_rsids, snps_by_rsid)
        six = timer()
        if not len(entries):
            raise NameError("Couldn not build covariance entries for (%s,%s)" %(name,gene))
        return entries,wat

    def updateFoundCovariance(self, gene, name):
        found = None
        if gene in self.found_genes_for_covariance:
            found = self.found_genes_for_covariance[gene]
            logging.info("Gene %s found again for %s", gene, name)
        else:
            found = []
            self.found_genes_for_covariance[gene] = found
        found.append(name)

    def buildRelatedData(self, rsids_from_genes, snps_by_rsid, weights_in_gene):
        related_rsids = []
        related_data = []
        l = len(rsids_from_genes)
        if self.max_snps_in_gene and l > self.max_snps_in_gene:
            logging.info("Skipping covariance too large: %d", l)
            return related_data, related_rsids

        for rsid in rsids_from_genes:
            if not rsid in snps_by_rsid:
                logging.log(5, "related rsid %s not present in genotype data", rsid)
                continue
            related_snp = snps_by_rsid[rsid]
            freq = sum(related_snp.data)*1.0/(2*len(related_snp.data))
            if self.min_maf_filter and self.min_maf_filter > freq:
                logging.log(6, "related rsid %s below min maf: %s", rsid, freq)
                continue

            if self.max_maf_filter and self.max_maf_filter < freq:
                logging.log(6, "related rsid %s  above max maf: %s", rsid, freq)
                continue
            related_rsids.append(rsid)
        def Getpos(elem):
            return snps_by_rsid[elem].position
        related_rsids_new = sorted(related_rsids, key = Getpos)
        for emen in related_rsids_new:
            related = snps_by_rsid[emen]
            data = related.data
            weight = weights_in_gene[emen]
            if weight.ref_allele == related.eff_allele and weight.eff_allele == related.ref_allele:
                logging.log(7, "related rsid %s has alleles flipped compared to model, transforming dosage", rsid)
                data = map(lambda x: 2-x, data)
            related_data.append(data)
        return related_rsids_new, related_data

    def buildMatrixOutputEntries(self, matrix, rsids_from_genes, related_rsids, snps_by_rsid):
        entries = []
        count  = 0
        for i in xrange(0, len(related_rsids)):
            rsid_i = related_rsids[i]
            value = str(matrix[i])
            entries.append((rsid_i, rsid_i, value))
        value = "NA"
        for k in xrange(0,len(related_rsids)):
            rsid_k = rsids_from_genes[k]
            if rsid_k not in related_rsids:
                entries.append((rsid_k,rsid_k,value))
        return entries

    def addToCorrelationFile(self, weight_db_logic, name, snps, snps_by_rsid):
        logging.info("Building correlation database for %s-%s", name, self.weight_db)
        genes = weight_db_logic.weights_by_gene.keys()
        total_genes = len(genes)
        last_reported_percent = 0
        processed = 0
        for gene in genes:
            processed += 1
            percent = int(processed*100.0 / total_genes)
            if percent == last_reported_percent+10:
                logging.info("%d percent genes processed", percent)
                last_reported_percent = percent

            entries = self.buildCorrelationEntries(name, gene, weight_db_logic, snps_by_rsid)

            if len(entries) == 0:
                logging.log(6,"Gene %s has no snps in current file", gene)
                continue

            self.addToFile(self.correlation_output, gene, entries)

    def buildCorrelationEntries(self, name, gene, weight_db_logic, snps_by_rsid):
        weights_in_gene = weight_db_logic.weights_by_gene[gene]
        rsids_from_genes = weights_in_gene.keys()

        #gather as much data as we can work on
        related_rsids, related_data = self.buildRelatedData(rsids_from_genes, snps_by_rsid, weights_in_gene)

        if len(related_rsids) == 0:
            return []

        self.updateFoundCorrelation(gene, name)

        #correlation matrix of related SNP's data
        array = numpy.array(related_data)
        cor = numpy.corrcoef(array)

        #translate into sql entries
        entries = self.buildMatrixOutputEntries(cor, rsids_from_genes, related_rsids, snps_by_rsid)
        if not len(entries):
            raise NameError("Couldn not build correlation entries for (%s,%s)" %(name,gene))
        return entries

    def updateFoundCorrelation(self, gene, name):
        found = None
        if gene in self.found_genes_for_correlation:
            found = self.found_genes_for_correlation[gene]
            logging.info("Gene %s found again for %s", gene, name)
        else:
            found = []
            self.found_genes_for_correlation[gene] = found
        found.append(name)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Build correlations and/or covariances from PHASE3 data and weights database.')

    parser.add_argument("--verbosity",
                        help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything",
                        default = "10")

    parser.add_argument("--weight_db",
                        help="name of weight db in data folder",
                        default="data/DGN-WB_0.5.db")

    parser.add_argument("--input_folder",
                        help="name of folder containing PHASE 3 data",
                        default="intermediate/TGF_EUR")

    parser.add_argument("--correlation_output",
                        help="Name of file to dump correlation results in.",
                        default=None)

    parser.add_argument("--covariance_output",
                        help="Name of file to dump covariance results in. Defaults to 'intermediate/cov/' + file name prefix from '--weight_db' argument",
                        default=None)

    parser.add_argument('--input_format',
                   help='Input dosage files format. Valid options are: IMPUTE, PrediXcan',
                   default=Formats.PrediXcan)

    parser.add_argument('--min_maf_filter',
                   help="Filter snps according to this maf",
                   default=None)

    parser.add_argument('--max_maf_filter',
                   help="Filter snps according to this maf",
                   default=None)

    parser.add_argument("--max_snps_in_gene",
                        help="Ignore any gene that has snps above this value",
                        type=int,
                        default=None)

    args = parser.parse_args()

    Logging.configureLogging(int(args.verbosity))

    work = ProcessWeightDB(args)
    work.run()
