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
public class TcgaSomaticMutationConverter extends BioDirectoryConverter
{
    //
    private static final String DATASET_TITLE = "TCGA Mutation Data";
    private static final String DATA_SOURCE_NAME = "TCGA Gene-level non-silent mutation PANCAN";

    private static final String TAXON_ID = "9606"; // Human Taxon ID

    private static final String MUTATION_TSV_FILE = "mc3.v0.2.8.PUBLIC.nonsilentGene.xena";

    private Map<String, String> genes = new HashMap<String, String>();
    private Map<String, String> samples = new HashMap<String, String>();

    protected IdResolver rslv;
    private static final Logger LOG = Logger.getLogger(TcgaSomaticMutationConverter.class);

    private String organismIdentifier; // Not the taxon ID. It references the object that is created into the database.
    /**
     * Constructor
     * @param writer the ItemWriter used to handle the resultant items
     * @param model the Model
     */
    public TcgaSomaticMutationConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
        if (rslv == null) {
            rslv = IdResolverService.getIdResolverByOrganism(TAXON_ID);
        }
    }

    /**
     * 
     *
     * {@inheritDoc}
     */
    public void process(File dataDir) throws Exception {
        Map<String, File> files = readFilesInDir(dataDir);

        organismIdentifier = getOrganism(TAXON_ID);

        processMutationData(new FileReader(files.get(MUTATION_TSV_FILE)));
    }

    private Map<String, File> readFilesInDir(File dir) {
        Map<String, File> files = new HashMap<String, File>();
        for (File file : dir.listFiles()) {
            files.put(file.getName(), file);
        }
        return files;
    }

    private void processMutationData(Reader reader) throws ObjectStoreException, IOException {
        Iterator<?> lineIter = FormattedTextParser.parseTabDelimitedReader(reader);
        // header
        String[] firstLine = (String[]) lineIter.next();
        ArrayList<String> samples = new ArrayList<String>();
        for(int i = 1; i < firstLine.length; i++) {
            String formattedSample = firstLine[i].split(" ")[0].trim();
            samples.add(formattedSample);
        }

        //lineIter.next();
        while (lineIter.hasNext()) {
            String[] line = (String[]) lineIter.next();

            String gene = line[0];
            for(int i = 1; i < line.length; i++) {
                String mutationValue = line[i];
                String theSampleForThisItem = samples.get(i-1);
                Item ExpressionItem;

                ExpressionItem = createItem("TCGAMutation");

                if(!gene.isEmpty()) {
                    String geneId = getGeneId(gene);

                    if (StringUtils.isEmpty(geneId)) {
                        continue;
                    }

                    ExpressionItem.setReference("gene", geneId);
                } else {
                    continue;
                }

                if(!theSampleForThisItem.isEmpty()) {
                    ExpressionItem.setReference("sampleID", getSample(theSampleForThisItem));
                } else {
                    continue;
                }

                if(!mutationValue.isEmpty() && StringUtils.isNumeric(mutationValue)) {
                    ExpressionItem.setAttribute("TcgaSomaticMutationValue", mutationValue);
                } else {
                    continue;
                }

                store(ExpressionItem);
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

    public String getSample(String identifier) {
        String refId = samples.get(identifier);
        if (refId == null) {
            Item cl = createItem("TCGASample");
            cl.setAttribute("SampleID", identifier);
            try {
                store(cl);
            } catch (ObjectStoreException e) {
                throw new RuntimeException("failed to store sample with ID: " + identifier, e);
            }
            refId = cl.getIdentifier();
            samples.put(identifier, refId);
        }
        return refId;
    }
}
