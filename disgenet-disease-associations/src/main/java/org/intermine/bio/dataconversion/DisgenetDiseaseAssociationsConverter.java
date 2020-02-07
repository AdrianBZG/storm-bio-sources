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

import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.util.FormattedTextParser;
import org.intermine.xml.full.Item;


/**
 * 
 * @author
 */
public class DisgenetDiseaseAssociationsConverter extends BioDirectoryConverter
{
    //
    private static final String DATASET_TITLE = "DisGeNET";
    private static final String DATA_SOURCE_NAME = "Curated gene-disease associations";

    private static final String TAXON_ID = "9606"; // Human Taxon ID

    private static final String DISGENET_FILE = "curated_gene_disease_associations.tsv";


    private Map<String, String> genes = new HashMap<String, String>();
    private Map<String, Item> diseases = new HashMap<String, Item>();


    private String organismIdentifier; // Not the taxon ID. It references the object that is created into the database.

    /**
     * Constructor
     * @param writer the ItemWriter used to handle the resultant items
     * @param model the Model
     */
    public DisgenetDiseaseAssociationsConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
    }

    @Override
    public void close() throws Exception {
        store(diseases.values());
    }

    /**
     * 
     *
     * {@inheritDoc}
     */
    public void process(File dataDir) throws Exception {

        Map<String, File> files = readFilesInDir(dataDir);

        organismIdentifier = getOrganism(TAXON_ID);

        processAssociations(new FileReader(files.get(DISGENET_FILE)));

    }

    private Map<String, File> readFilesInDir(File dir) {
        Map<String, File> files = new HashMap<String, File>();
        for (File file : dir.listFiles()) {
            files.put(file.getName(), file);
        }
        return files;
    }

    private void processAssociations(Reader reader) throws ObjectStoreException, IOException {
        Iterator<?> lineIter = FormattedTextParser.parseTabDelimitedReader(reader);
        // Skip header
        lineIter.next();

        while (lineIter.hasNext()) {
            String[] line = (String[]) lineIter.next();

            String geneSymbol = line[0];
            String diseaseId = line[4];
            String diseaseName = line[5];
            String diseaseType = line[6];
            String score = line[9];

            Item disease = getDisease(diseaseId);

            if (disease == null) {
                disease = createItem("Disease");
                disease.setAttribute("diseaseId", diseaseId);
                disease.setAttribute("diseaseName", diseaseName);
                disease.setAttribute("diseaseType", diseaseType);
                store(disease);
            }

            Item interactionItem;
            interactionItem = createItem("DiseaseAssociation");
            interactionItem.setReference("gene", getGene(geneSymbol));
            interactionItem.setReference("disease", disease);
            interactionItem.setAttribute("associationScore", score);
            store(interactionItem);
        }
    }

    public String getGene(String identifier) {
        String refId = genes.get(identifier);
        if (refId == null) {
            Item gene = createItem("Gene");
            gene.setAttribute("symbol", identifier);
            gene.setReference("organism", getOrganism(TAXON_ID));
            try {
                store(gene);
            } catch (ObjectStoreException e) {
                throw new RuntimeException("failed to store gene with primary identifier: " + identifier, e);
            }
            refId = gene.getIdentifier();
            genes.put(identifier, refId);
        }
        return refId;
    }

    private Item getDisease(String diseaseId) {
        Item disease = diseases.get(diseaseId);
        return disease;
    }
}
