"""
CS123A project that takes an inputted sequence and returns the index of the beginning of the predicted promoter site
"""

def promoter_prediction(text_file):

#reads sequence text file, skips first line, and stores in a one line string
    with open(text_file) as f:
        next(f)
        seq = f.read().replace('\n', '')

#list containing all possible Inr sequences (consensus seq = BBCABW)
    Inr_list = ["CCCACA", "CCCACT", "CCCAGA", "CCCAGT", "CCCATA", "CCCATT", "CGCACA", "CGCACT", "CGCAGA", "CGCAGT",
           "CGCATA","CGCATT", "CTCACA", "CTCACT", "CTCAGA", "CTCAGT", "CTCATA", "CTCATT", "GCCACA", "GCCACT", "GCCAGA",
            "GCCAGT", "GCCATA", "GCCATT", "GGCACA", "GGCACT", "GGCAGA", "GGCAGT", "GGCATA", "GGCATT", "GTCACA", "GTCACT",
            "GTCAGA", "GTCAGT", "GTCATA", "GTCATT", "TCCACA", "TCCACT", "TCCAGA", "TCCAGT", "TCCATA", "TCCATT", "TGCACA",
            "TGCACT", "TGCAGA", "TGCAGT", "TGCATA", "TGCATT", "TTCACA", "TTCACT", "TTCAGA", "TTCAGT", "TTCATA", "TTCATT"]

#list to store predicted promoter start sites
    pred_prom = []

#dictionary to store predicted promoter start sites and the number of G and C characters counted
    pred_prom_gc_count = {}

#searches for possible Inr sequences in inputted sequence and returns index
    for subseq in Inr_list:

        #if subsequence is not found, then that element is skipped
        if seq.find(subseq) == -1:
            pass
        #if Inr element is located, it then looks for TATA box within 30 bp of the Inr element
        else:
            #store index of where the inr sequence begins in the inputted sequence
            pot_promoter = seq.find(subseq)

            #looking for the TATA box sequence within 30 bp of the inr sequence
            pot_prom_area = seq[(seq.find(subseq) - 30):(seq.find(subseq) + 30)]

            #if the TATA box is found, then the index of the predicted promoter sequence is added to predicted promoter list
            if pot_prom_area.find("TATA") != -1:
                pred_prom.append(pot_promoter)

            #counts the number of CGs that occur in predicted promoter sites that are close to a TATA box
                gc_count = pot_prom_area.count("CG")
                pred_prom_gc_count[pot_promoter] = gc_count

    #sort list of predicted promoter indexes in ascending order
    pred_prom.sort()

    #find predicted promoter index that had highest number of GC sequences, which could be indicative of CpG islands
    best_prom = max(pred_prom_gc_count, key=pred_prom_gc_count.get)

    #prints all the potential promoter sites, and the single best predicted promoter site
    print("")
    print("For the gene", text_file)
    print("The predicted promoter start sites are:", ", ".join(map(str, pred_prom)))
    print("The best predicted promoter site is:", best_prom)

promoter_prediction("HBB.txt")
promoter_prediction("COVID-19.txt")
promoter_prediction("CFTR.txt")
promoter_prediction("TNF.txt")
promoter_prediction("BRCA1.txt")
promoter_prediction("PKD2.txt")





