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
public class TcgaSampleMetadataConverter extends BioDirectoryConverter
{
    //
    private static final String DATASET_TITLE = "TCGA Sample Metadata";
    private static final String DATA_SOURCE_NAME = "TCGA sample type and primary disease (PANCAN)";

    private static final String TAXON_ID = "9606"; // Human Taxon ID

    private static final String SAMPLE_INFO_CSV_FILE = "TCGA_phenotype_denseDataOnlyDownload.tsv";

    private Map<String, String> samples = new HashMap<String, String>();

    private String organismIdentifier; // Not the taxon ID. It references the object that is created into the database.

    /**
     * Constructor
     * @param writer the ItemWriter used to handle the resultant items
     * @param model the Model
     */
    public TcgaSampleMetadataConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
    }

    /**
     * 
     *
     * {@inheritDoc}
     */
    public void process(File dataDir) throws Exception {

        Map<String, File> files = readFilesInDir(dataDir);

        organismIdentifier = getOrganism(TAXON_ID);

        processSamples(new FileReader(files.get(SAMPLE_INFO_CSV_FILE)));

    }

    private Map<String, File> readFilesInDir(File dir) {
        Map<String, File> files = new HashMap<String, File>();
        for (File file : dir.listFiles()) {
            files.put(file.getName(), file);
        }
        return files;
    }

    private void processSamples(Reader reader) throws ObjectStoreException, IOException {
        Iterator<?> lineIter = FormattedTextParser.parseTabDelimitedReader(reader);
        // Skip header
        lineIter.next();
        while (lineIter.hasNext()) {
            String[] line = (String[]) lineIter.next();

            String SampleID = line[0];

            if(samples.containsKey(SampleID)) {
                continue;
            }

            String SampleTypeID = line[1];
            String SampleType = line[2];
            String PrimaryDisease = line[3];

            Item sampleItem;

            sampleItem = createItem("TCGASample");

            if(!SampleID.isEmpty()) {
                sampleItem.setAttribute("SampleID", SampleID);
            } else {
                continue;
            }

            if(!SampleTypeID.isEmpty()) {
                sampleItem.setAttribute("SampleTypeID", SampleTypeID);
            } else {
                sampleItem.setAttribute("SampleTypeID", "Not specified");
            }

            if(!SampleType.isEmpty()) {
                sampleItem.setAttribute("SampleType", SampleType);
            } else {
                sampleItem.setAttribute("SampleType", "Not specified");
            }

            if(!PrimaryDisease.isEmpty()) {
                sampleItem.setAttribute("PrimaryDisease", PrimaryDisease);
            } else {
                sampleItem.setAttribute("PrimaryDisease", "Not specified");
            }

            store(sampleItem);
            samples.put(SampleID, sampleItem.getIdentifier());
        }
    }
}
