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

import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.util.FormattedTextParser;
import org.intermine.xml.full.Item;


/**
 * 
 * @author
 */
public class DepmapCcleMutationsConverter extends BioDirectoryConverter
{
    //
    private static final String DATASET_TITLE = "DepMap CCLE Mutations Data";
    private static final String DATA_SOURCE_NAME = "DepMap Public 19Q3";

    private static final String TAXON_ID = "9606"; // Human Taxon ID

    private static final String MUTATIONS_CSV_FILE = "CCLE_mutations.csv";

    private Map<String, String> genes = new HashMap<String, String>();
    private Map<String, String> cellLines = new HashMap<String, String>();

    private String organismIdentifier; // Not the taxon ID. It references the object that is created into the database.

    /**
     * Constructor
     * @param writer the ItemWriter used to handle the resultant items
     * @param model the Model
     */
    public DepmapCcleMutationsConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
    }

    public void process(File dataDir) throws Exception {
        Map<String, File> files = readFilesInDir(dataDir);

        organismIdentifier = getOrganism(TAXON_ID);

        processMutationsData(new FileReader(files.get(MUTATIONS_CSV_FILE)));
    }

    private Map<String, File> readFilesInDir(File dir) {
        Map<String, File> files = new HashMap<String, File>();
        for (File file : dir.listFiles()) {
            files.put(file.getName(), file);
        }
        return files;
    }

    private void processMutationsData(Reader reader) throws ObjectStoreException, IOException {
        Iterator<?> lineIter = FormattedTextParser.parseCsvDelimitedReader(reader);

        // Skip header
        lineIter.next();

        while (lineIter.hasNext()) {
            String[] line = (String[]) lineIter.next();

            Item mutationsItem;

            mutationsItem = createItem("DepMapMutations");

            /*
                -	Hugo Symbol 0
                -	Chromosome 3
                -	Start position 4
                -	End position 5
                -	Strand 6
                -	Variant classification 7
                -	Variant type 8
                -	Genome change 13
                -	Annotation transcript 14
                -	isDeleterious 19
                -	isTCGAhotspot 20
                -	TCGAhsCnt 21
                -	isCOSMIChotspot 22
                -	COSMIChsCnt 23
                -	Variant annotation 32
                -	DepMap ID 33
            */

            String hugoSymbol = line[0];
            String chromosome = line[3];
            String startPosition = line[4];
            String endPosition = line[5];
            String strand = line[6];
            String variantClassification = line[7];
            String variantType = line[8];
            String genomeChange = line[13];
            String annotationTranscript = line[14];
            String isDeleterious = line[19];
            String isTCGAhotspot = line[20];
            String TCGAhsCnt = line[21];
            String isCOSMIChotspot = line[22];
            String COSMIChsCnt = line[23];
            String variantAnnotation = line[32];
            String DepMapID = line[33];



            if(!DepMapID.isEmpty()) {
                mutationsItem.setReference("depMapID", getCellLine(DepMapID));
            } else {
                continue;
            }

            if(!hugoSymbol.isEmpty()) {
                mutationsItem.setReference("gene", getGene(hugoSymbol));
            } else {
                continue;
            }

            if(!chromosome.isEmpty()) {
                mutationsItem.setAttribute("Chromosome", chromosome);
            } else {
                mutationsItem.setAttribute("Chromosome", "Not specified");
            }

            if(!startPosition.isEmpty()) {
                mutationsItem.setAttribute("Start", startPosition);
            } else {
                mutationsItem.setAttribute("Start", "Not specified");
            }

            if(!endPosition.isEmpty()) {
                mutationsItem.setAttribute("End", endPosition);
            } else {
                mutationsItem.setAttribute("End", "Not specified");
            }

            if(!strand.isEmpty()) {
                mutationsItem.setAttribute("Strand", strand);
            } else {
                mutationsItem.setAttribute("Strand", "Not specified");
            }

            if(!variantClassification.isEmpty()) {
                mutationsItem.setAttribute("VariantClassification", variantClassification);
            } else {
                mutationsItem.setAttribute("VariantClassification", "Not specified");
            }

            if(!variantType.isEmpty()) {
                mutationsItem.setAttribute("VariantType", variantType);
            } else {
                mutationsItem.setAttribute("VariantType", "Not specified");
            }

            if(!genomeChange.isEmpty()) {
                mutationsItem.setAttribute("GenomeChange", genomeChange);
            } else {
                mutationsItem.setAttribute("GenomeChange", "Not specified");
            }

            if(!annotationTranscript.isEmpty()) {
                mutationsItem.setAttribute("AnnotationTranscript", annotationTranscript);
            } else {
                mutationsItem.setAttribute("AnnotationTranscript", "Not specified");
            }

            if(!isDeleterious.isEmpty()) {
                mutationsItem.setAttribute("isDeleterious", isDeleterious);
            } else {
                mutationsItem.setAttribute("isDeleterious", "Not specified");
            }

            if(!isTCGAhotspot.isEmpty()) {
                mutationsItem.setAttribute("isTCGAhotspot", isTCGAhotspot);
            } else {
                mutationsItem.setAttribute("isTCGAhotspot", "Not specified");
            }

            if(!TCGAhsCnt.isEmpty()) {
                mutationsItem.setAttribute("TCGAhsCnt", TCGAhsCnt);
            } else {
                mutationsItem.setAttribute("TCGAhsCnt", "Not specified");
            }

            if(!isCOSMIChotspot.isEmpty()) {
                mutationsItem.setAttribute("isCOSMIChotspot", isCOSMIChotspot);
            } else {
                mutationsItem.setAttribute("isCOSMIChotspot", "Not specified");
            }

            if(!COSMIChsCnt.isEmpty()) {
                mutationsItem.setAttribute("COSMIChsCnt", COSMIChsCnt);
            } else {
                mutationsItem.setAttribute("COSMIChsCnt", "Not specified");
            }

            if(!variantAnnotation.isEmpty()) {
                mutationsItem.setAttribute("VariantAnnotation", variantAnnotation);
            } else {
                mutationsItem.setAttribute("VariantAnnotation", "Not specified");
            }

            store(mutationsItem);
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

    public String getCellLine(String identifier) {
        String refId = cellLines.get(identifier);
        if (refId == null) {
            Item cl = createItem("CellLine");
            cl.setAttribute("DepMapID", identifier);
            try {
                store(cl);
            } catch (ObjectStoreException e) {
                throw new RuntimeException("failed to store cell line with DepMapID: " + identifier, e);
            }
            refId = cl.getIdentifier();
            cellLines.put(identifier, refId);
        }
        return refId;
    }
}
