package org.intermine.bio.dataconversion;

/*
 * Copyright (C) 2002-2018 FlyMine
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  See the LICENSE file for more
 * information or http://www.gnu.org/copyleft/lesser.html.
 *
 */

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.util.FormattedTextParser;
import org.intermine.xml.full.Item;


/**
 * 
 * @author
 */
public class StormTargetsAnalysesConverter extends BioDirectoryConverter
{
    //
    private static final String DATASET_TITLE = "STORM Targets Analyses";
    private static final String DATA_SOURCE_NAME = "Results for the 2020 target triaging analyses on STORM Targets";

    private static final String TAXON_ID = "9606"; // Human Taxon ID

    private static final String TARGETS_ANALYSES_FILE = "DepMap_RME_results_with_outliers.csv";

    private Map<String, String> genes = new HashMap<String, String>();

    protected IdResolver rslv;
    private static final Logger LOG = Logger.getLogger(StormTargetsAnalysesConverter.class);

    private String organismIdentifier;
    public StormTargetsAnalysesConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
        if (rslv == null) {
            rslv = IdResolverService.getIdResolverByOrganism(TAXON_ID);
        }
    }

    private boolean isDouble(String str) {
        try {
            // check if it can be parsed as any double
            double x = Double.parseDouble(str);
            // check if the double can be converted without loss to an int
            if (x == (int) x)
                // if yes, this is an int, thus return false
                return false;
            // otherwise, this cannot be converted to an int (e.g. "1.2")
            return true;
            // short version: return x != (int) x;
        }
        catch(NumberFormatException e) {
            return false;
        }
    }

    public void process(File dataDir) throws Exception {
        Map<String, File> files = readFilesInDir(dataDir);

        organismIdentifier = getOrganism(TAXON_ID);

        processTargetsAnalyses(new FileReader(files.get(TARGETS_ANALYSES_FILE)));
    }

    private Map<String, File> readFilesInDir(File dir) {
        Map<String, File> files = new HashMap<String, File>();
        for (File file : dir.listFiles()) {
            files.put(file.getName(), file);
        }
        return files;
    }

    private void processTargetsAnalyses(Reader reader) throws ObjectStoreException, IOException {
        Iterator<?> lineIter = FormattedTextParser.parseCsvDelimitedReader(reader);
        // header
        String[] firstLine = (String[]) lineIter.next();

        //lineIter.next();
        while (lineIter.hasNext()) {
            String[] line = (String[]) lineIter.next();

            String gene = line[0];
            String broadMedian = line[6];
            String broadEffectFraction = line[7];
            String broadCommonEssential = line[8];
            String broadSkewedLrt = line[9];
            String sangerMedian = line[10];
            String sangerEffectFraction = line[11];
            String sangerCommonEssential = line[12];
            String sangerSkewedLrt = line[13];
            String shrnaMedian = line[14];
            String shrnaEffectFraction = line[15];
            String shrnaCommonEssential = line[16];
            String shrnaSkewedLrt = line[17];
            String broadSangerCor = line[18];
            String broadSangerDiff = line[19];
            String broadShrnaCor = line[20];
            String broadShrnaDiff = line[21];
            String sangerShrnaCor = line[22];
            String sangerShrnaDiff = line[23];
            String broadOutliersCount = line[24];
            String broadOutliersMeanZscore = line[25];
            String broadOutliersCellLines = line[26];
            String broadOutliersTopLineage = line[27];
            String broadOutliersLineageCount = line[28];
            String broadOutliersLineagePvalue = line[29];
            String broadOutliersLineageQvalue = line[30];
            String sangerOutliersCount = line[31];
            String sangerOutliersMeanZscore = line[32];
            String sangerOutliersCellLines = line[33];
            String sangerOutliersTopLineage = line[34];
            String sangerOutliersLineageCount = line[35];
            String sangerOutliersLineagePvalue = line[36];
            String sangerOutliersLineageQvalue = line[37];
            String shrnaOutliersCount = line[38];
            String shrnaOutliersMeanZscore = line[39];
            String shrnaOutliersCellLines = line[40];
            String shrnaOutliersTopLineage = line[41];
            String shrnaOutliersLineageCount = line[42];
            String shrnaOutliersLineagePvalue = line[43];
            String shrnaOutliersLineageQvalue = line[44];


            Item IntegratedItem;


            String geneId = "";

            if(!gene.isEmpty()) {
                geneId = getGeneId(gene);

                if (StringUtils.isEmpty(geneId)) {
                    continue;
                }
            } else {
                continue;
            }

            if(!StringUtils.isEmpty(broadMedian) && isDouble(broadMedian)) {
                IntegratedItem = createItem("STORMTargetAnalyses");
                IntegratedItem.setReference("gene", geneId);
                IntegratedItem.setAttribute("screen", "broad");
                IntegratedItem.setAttribute("median", broadMedian);
                if(isDouble(sangerEffectFraction)) IntegratedItem.setAttribute("effectFraction", broadEffectFraction);
                if(isDouble(broadCommonEssential)) IntegratedItem.setAttribute("commonEssential", broadCommonEssential);
                if(isDouble(broadSkewedLrt)) IntegratedItem.setAttribute("skewedLrt", broadSkewedLrt);
                IntegratedItem.setAttribute("outliersCount", broadOutliersCount);
                if(isDouble(broadOutliersMeanZscore)) IntegratedItem.setAttribute("outliersMeanZscore", broadOutliersMeanZscore);
                if(!StringUtils.isEmpty(broadOutliersCellLines)) IntegratedItem.setAttribute("outliersCellLines", broadOutliersCellLines);
                if(!StringUtils.isEmpty(broadOutliersTopLineage)) IntegratedItem.setAttribute("outliersTopLineage", broadOutliersTopLineage);
                IntegratedItem.setAttribute("outliersLineageCount", broadOutliersLineageCount);
                if(isDouble(broadOutliersLineagePvalue)) IntegratedItem.setAttribute("outliersLineagePvalue", broadOutliersLineagePvalue);
                if(isDouble(broadOutliersLineageQvalue)) IntegratedItem.setAttribute("outliersLineageQvalue", broadOutliersLineageQvalue);
                if(isDouble(broadSangerCor)) IntegratedItem.setAttribute("broadSangerCor", broadSangerCor);
                if(isDouble(broadSangerDiff)) IntegratedItem.setAttribute("broadSangerDiff", broadSangerDiff);
                if(isDouble(broadShrnaCor)) IntegratedItem.setAttribute("broadShrnaCor", broadShrnaCor);
                if(isDouble(broadShrnaDiff)) IntegratedItem.setAttribute("broadShrnaDiff", broadShrnaDiff);
                if(isDouble(sangerShrnaCor)) IntegratedItem.setAttribute("sangerShrnaCor", sangerShrnaCor);
                if(isDouble(sangerShrnaDiff)) IntegratedItem.setAttribute("sangerShrnaDiff", sangerShrnaDiff);
                store(IntegratedItem);
            }

            if(!StringUtils.isEmpty(sangerMedian) && isDouble(sangerMedian)) {
                IntegratedItem = createItem("STORMTargetAnalyses");
                IntegratedItem.setReference("gene", geneId);
                IntegratedItem.setAttribute("screen", "sanger");
                IntegratedItem.setAttribute("median", sangerMedian);
                if(isDouble(sangerEffectFraction)) IntegratedItem.setAttribute("effectFraction", sangerEffectFraction);
                if(isDouble(sangerCommonEssential)) IntegratedItem.setAttribute("commonEssential", sangerCommonEssential);
                if(isDouble(sangerSkewedLrt)) IntegratedItem.setAttribute("skewedLrt", sangerSkewedLrt);
                IntegratedItem.setAttribute("outliersCount", sangerOutliersCount);
                if(isDouble(sangerOutliersMeanZscore)) IntegratedItem.setAttribute("outliersMeanZscore", sangerOutliersMeanZscore);
                if(!StringUtils.isEmpty(sangerOutliersCellLines)) IntegratedItem.setAttribute("outliersCellLines", sangerOutliersCellLines);
                if(!StringUtils.isEmpty(sangerOutliersTopLineage)) IntegratedItem.setAttribute("outliersTopLineage", sangerOutliersTopLineage);
                IntegratedItem.setAttribute("outliersLineageCount", sangerOutliersLineageCount);
                if(isDouble(sangerOutliersLineagePvalue)) IntegratedItem.setAttribute("outliersLineagePvalue", sangerOutliersLineagePvalue);
                if(isDouble(sangerOutliersLineageQvalue)) IntegratedItem.setAttribute("outliersLineageQvalue", sangerOutliersLineageQvalue);
                if(isDouble(broadSangerCor)) IntegratedItem.setAttribute("broadSangerCor", broadSangerCor);
                if(isDouble(broadSangerDiff)) IntegratedItem.setAttribute("broadSangerDiff", broadSangerDiff);
                if(isDouble(broadShrnaCor)) IntegratedItem.setAttribute("broadShrnaCor", broadShrnaCor);
                if(isDouble(broadShrnaDiff)) IntegratedItem.setAttribute("broadShrnaDiff", broadShrnaDiff);
                if(isDouble(sangerShrnaCor)) IntegratedItem.setAttribute("sangerShrnaCor", sangerShrnaCor);
                if(isDouble(sangerShrnaDiff)) IntegratedItem.setAttribute("sangerShrnaDiff", sangerShrnaDiff);
                store(IntegratedItem);
            }

            if(!StringUtils.isEmpty(shrnaMedian) && isDouble(shrnaMedian)) {
                IntegratedItem = createItem("STORMTargetAnalyses");
                IntegratedItem.setReference("gene", geneId);
                IntegratedItem.setAttribute("screen", "shrna");
                IntegratedItem.setAttribute("median", shrnaMedian);
                if(isDouble(shrnaEffectFraction)) IntegratedItem.setAttribute("effectFraction", shrnaEffectFraction);
                if(isDouble(shrnaCommonEssential)) IntegratedItem.setAttribute("commonEssential", shrnaCommonEssential);
                if(isDouble(shrnaSkewedLrt)) IntegratedItem.setAttribute("skewedLrt", shrnaSkewedLrt);
                IntegratedItem.setAttribute("outliersCount", shrnaOutliersCount);
                if(isDouble(shrnaOutliersMeanZscore)) IntegratedItem.setAttribute("outliersMeanZscore", shrnaOutliersMeanZscore);
                if(!StringUtils.isEmpty(shrnaOutliersCellLines)) IntegratedItem.setAttribute("outliersCellLines", shrnaOutliersCellLines);
                if(!StringUtils.isEmpty(shrnaOutliersTopLineage)) IntegratedItem.setAttribute("outliersTopLineage", shrnaOutliersTopLineage);
                IntegratedItem.setAttribute("outliersLineageCount", shrnaOutliersLineageCount);
                if(isDouble(shrnaOutliersLineagePvalue)) IntegratedItem.setAttribute("outliersLineagePvalue", shrnaOutliersLineagePvalue);
                if(isDouble(shrnaOutliersLineageQvalue)) IntegratedItem.setAttribute("outliersLineageQvalue", shrnaOutliersLineageQvalue);
                if(isDouble(broadSangerCor)) IntegratedItem.setAttribute("broadSangerCor", broadSangerCor);
                if(isDouble(broadSangerDiff)) IntegratedItem.setAttribute("broadSangerDiff", broadSangerDiff);
                if(isDouble(broadShrnaCor)) IntegratedItem.setAttribute("broadShrnaCor", broadShrnaCor);
                if(isDouble(broadShrnaDiff)) IntegratedItem.setAttribute("broadShrnaDiff", broadShrnaDiff);
                if(isDouble(sangerShrnaCor)) IntegratedItem.setAttribute("sangerShrnaCor", sangerShrnaCor);
                if(isDouble(sangerShrnaDiff)) IntegratedItem.setAttribute("sangerShrnaDiff", sangerShrnaDiff);
                store(IntegratedItem);
            }


        }
    }

    private String getGeneId(String primaryIdentifier) throws ObjectStoreException {
        String resolvedIdentifier = resolveGene(primaryIdentifier);
        if (StringUtils.isEmpty(resolvedIdentifier)) {
            return null;
        }
        String geneId = genes.get(resolvedIdentifier);
        if (geneId == null) {
            Item gene = createItem("Gene");
            gene.setAttribute("primaryIdentifier", resolvedIdentifier);
            //gene.setAttribute("symbol", primaryIdentifier);
            //gene.setReference("organism", getOrganism(TAXON_ID));
            store(gene);
            geneId = gene.getIdentifier();
            genes.put(resolvedIdentifier, geneId);
        }
        return geneId;
    }

    private String resolveGene(String identifier) {
        String id = identifier;

        if (rslv != null && rslv.hasTaxon(TAXON_ID)) {
            int resCount = rslv.countResolutions(TAXON_ID, identifier);
            if (resCount != 1) {
                LOG.info("RESOLVER: failed to resolve gene to one identifier, ignoring gene: "
                        + identifier + " count: " + resCount + " Human identifier: "
                        + rslv.resolveId(TAXON_ID, identifier));
                return null;
            }
            id = rslv.resolveId(TAXON_ID, identifier).iterator().next();
        }
        return id;
    }
}
