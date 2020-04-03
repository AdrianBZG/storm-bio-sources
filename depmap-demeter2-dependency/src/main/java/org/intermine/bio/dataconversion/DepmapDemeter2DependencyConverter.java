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
import java.util.ArrayList;
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
public class DepmapDemeter2DependencyConverter extends BioDirectoryConverter
{
    //
    private static final String DATASET_TITLE = "DepMap DEMETER2 Gene Dependency";
    private static final String DATA_SOURCE_NAME = "DepMap DEMETER2 Gene Dependency";

    private static final String TAXON_ID = "9606"; // Human Taxon ID

    private static final String CN_CSV_FILE = "D2_combined_gene_dep_scores.csv";

    protected IdResolver rslv;
    private static final Logger LOG = Logger.getLogger(DepmapDemeter2DependencyConverter.class);

    private Map<String, String> genes = new HashMap<String, String>();
    private Map<String, String> cellLinesMap = new HashMap<String, String>();

    private String organismIdentifier; // Not the taxon ID. It references the object that is created into the database.


    /**
     * Constructor
     * @param writer the ItemWriter used to handle the resultant items
     * @param model the Model
     */
    public DepmapDemeter2DependencyConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
        if (rslv == null) {
            rslv = IdResolverService.getIdResolverByOrganism(TAXON_ID);
        }
    }

    public void process(File dataDir) throws Exception {

        Map<String, File> files = readFilesInDir(dataDir);

        organismIdentifier = getOrganism(TAXON_ID);

        processDependency(new FileReader(files.get(CN_CSV_FILE)));

    }

    private Map<String, File> readFilesInDir(File dir) {
        Map<String, File> files = new HashMap<String, File>();
        for (File file : dir.listFiles()) {
            files.put(file.getName(), file);
        }
        return files;
    }

    private void processDependency(Reader reader) throws ObjectStoreException, IOException {
        Iterator<?> lineIter = FormattedTextParser.parseCsvDelimitedReader(reader);
        // header
        String[] firstLine = (String[]) lineIter.next();
        ArrayList<String> cellLines = new ArrayList<String>();
        for(int i = 1; i < firstLine.length; i++) {
            String formattedCL = firstLine[i].trim();
            cellLines.add(formattedCL);
        }

        //lineIter.next();
        while (lineIter.hasNext()) {
            String[] line = (String[]) lineIter.next();

            String gene = line[0].split(" ")[0].trim();
            for(int i = 1; i < line.length; i++) {
                String dependencyValue = line[i];
                String theCLForThisItem = cellLines.get(i-1);
                Item DEMETER2Item;

                DEMETER2Item = createItem("DEMETER2Dependency");

                if(!gene.isEmpty()) {
                    DEMETER2Item.setReference("CCLEName", getCellLine(theCLForThisItem));
                } else {
                    continue;
                }

                if(!gene.isEmpty()) {
                    String geneId = getGeneId(gene);

                    if (StringUtils.isEmpty(geneId)) {
                        continue;
                    }

                    DEMETER2Item.setReference("gene", geneId);
                } else {
                    continue;
                }

                if(!dependencyValue.isEmpty() && StringUtils.isNumeric(dependencyValue)) {
                    DEMETER2Item.setAttribute("DepMapDEMETER2DependencyValue", dependencyValue);
                } else {
                    continue;
                }

                store(DEMETER2Item);
                //cellLines.put(cellLine, CopyNumberItem.getIdentifier());
            }
        }
    }

    private String getGeneId(String symbol) throws ObjectStoreException {
        //String resolvedIdentifier = resolveGene(primaryIdentifier);
        /*if (StringUtils.isEmpty(resolvedIdentifier)) {
            return null;
        }*/
        String geneId = genes.get(symbol);
        if (geneId == null) {
            Item gene = createItem("Gene");
            gene.setAttribute("symbol", symbol);
            gene.setReference("organism", getOrganism(TAXON_ID));
            store(gene);
            geneId = gene.getIdentifier();
            genes.put(symbol, geneId);
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

    public String getCellLine(String identifier) {
        String refId = cellLinesMap.get(identifier);
        if (refId == null) {
            Item cl = createItem("CellLine");
            cl.setAttribute("CCLEname", identifier);
            try {
                store(cl);
            } catch (ObjectStoreException e) {
                throw new RuntimeException("failed to store cell line with CCLE name: " + identifier, e);
            }
            refId = cl.getIdentifier();
            cellLinesMap.put(identifier, refId);
        }
        return refId;
    }
}
