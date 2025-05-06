import requests
import time
import re
from Bio import Entrez, SeqIO
from xml.etree import ElementTree as ET
from collections import Counter
import os
from urllib.request import urlretrieve
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.chrome.service import Service
import time
import csv
import pandas as pd
from collections import defaultdict
from pathlib import Path
import google.generativeai as genai
from google.api_core.exceptions import ResourceExhausted

BASE_DIR = os.getcwd()

api_key = os.getenv("NCBI_API_KEY")


input_seq = input("Enter the 4-letter PDB ID: ").strip()
URL = f"https://www.rcsb.org/fasta/entry/{input_seq}"
response = requests.get(URL)
print(response.text)


if response.status_code == 200:
    fasta_lines = response.text.strip().split("\n")
    chain_groups = {}  
    current_chain_group = None
    sequence = ""
    
    
    for line in fasta_lines:
        if line.startswith(">"): 
            if current_chain_group and sequence:
                chain_groups[current_chain_group]["sequence"] = sequence #you take the current chain sequence and append it to the empty sequence of the first entry
            
            parts = line.split("|")
            #print(parts)
            if len(parts) > 2:
                chain_group_id = parts[0][1:].strip() #parts[0]--> goes to the chain say ">1F1B_1" and in that it takes from 1 till the end(basically it takes each letter and does slicing of the string) 
                chains = parts[1].replace("Chains ", "").split(", ")  # replacing the word Chains with nothing and then splitting based on 
                title = parts[2].strip()  # Protein title
                organism_info = parts[3].strip()
                #print(organism_info.split())

                # Extract organism name and taxonomic ID
                temp=organism_info.split()
                #print(temp)
                temp_1=temp[:-1]
                org_name=" ".join(temp_1)
                
                
                #taxid extraction
                temp_1=temp[-1]
                tax_id=temp_1.lstrip("(").rstrip(")")
                  
                # Store in dictionary
                chain_groups[chain_group_id] = { # A Dictionary in a dictionary
                    "chains": chains,
                    "title": title,
                    "organism": org_name,
                    "taxid": tax_id,
                    "sequence": ""
                }
                current_chain_group = chain_group_id  # Set current group
            sequence = ""  
        else:
            sequence += line.strip()  
    
    if current_chain_group and sequence:
        chain_groups[current_chain_group]["sequence"] = sequence # you again do the same for the last entry in the fasta file, if any, cuz the loop doesn't go back anymore
    #once you are finally out of the loop you see if there are multiple chains by simply counting the number of dictionaries(since it is a dict in a dict)
    if len(chain_groups) > 1:
        print("\nThere are multiple chains available for this PDB entry.")
        print("Please select one of the available chain groups:")
        #print(chain_groups.values())--> dict_keys(['1F1B_1', '1F1B_2'])
        #print(chain_groups.keys())-->dict_values([{'chains': ['A', 'C'], 'title': 'ASPARTATE CARBAMOYLTRANSFERASE CATALYTIC CHAIN', 'organism': 'Escherichia coli', 'taxid': '562', 'sequence': 'ANPLYQKHIISINDLSRDDLNLVLATAAKLKANPQPELLKHKVIASCFFEASTRTRLSFETSMHRLGASVVGFSDSANTSLGKKGETLADTISVISTYVDAIVMRHPQEGAARLATEFSGNVPVLNAGDGSNQHPTQTLLDLFTIQETQGRLDNLHVAMVGDLKYGRTVHSLTQALAKFDGNRFYFIAPDALAMPQYILDMLDEKGIAWSLHSSIEEVMAEVDILYMTRVQKERLDPSEYANVKAQFVLRASDLHNAKANMKVLHPLARVDEIATDVDKTPHAWYFQQAGNGIFARQALLALVLNRDLVL'}, {'chains': ['B', 'D'], 'title': 'ASPARTATE CARBAMOYLTRANSFERASE REGULATORY CHAIN', 'organism': 'Escherichia coli', 'taxid': '562', 'sequence': 'MTHDNKLQVEAIKRGTVIDHIPAQIGFKLLSLFKLTETDQRITIGLNLPSGEMGRKDLIKIENTFLSEDQVDQLALYAPQATVNRIDNYEVVGKSRPSLPERIDNVLVCPNSNCISHAEPVSSSFAVRKRANDIALKCKYCEKEFSHNVVLAN'}])
        #keys are the chain_ids like 1F1B_1 and 1F1B_2 
        #values are the dictionary which again contains keys and values
        # Display 
        index = 1
        for chain_group in chain_groups:#chain group refers to one dictionary in then entire dict-chain_groups
            chains = ", ".join(chain_groups[chain_group]["chains"])#you go to the main dictionary and then the sub dictionary and then the key "chains" and get its corresponding value
            print(f"{index}. {chain_group} (Chains: {chains})") # output-->1. 1F1B_1 (Chains: A, C)
            index += 1 #then increase the index so that you can access the second chain_id


        while True:
            try:
                choice = int(input("\nEnter the number corresponding to your desired chain group: "))
                if 1 <= choice <= len(chain_groups):
                    selected_chain_group = list(chain_groups.keys())[choice - 1]#choice-1 cuz of zero based indexing of list slicing, basically you are getting the chain_id
                    break
                else:
                    print("Invalid choice. Please select a valid number.")
            except ValueError:
                print("Invalid input. Please enter a number.")

    else:
        # Only one chain, select it automatically
        selected_chain_group = list(chain_groups.keys())[0]

    # Extract selected chain's information
    selected_info = chain_groups[selected_chain_group]

    #save taxonomy and organism name and it's genus to a separate variable
    organism_for_tblastn=selected_info["organism"]
    genus=organism_for_tblastn.split()[0]
    taxid=selected_info['taxid']
    title=selected_info['title']
    sequence = selected_info["sequence"]
    pdb_length=len(sequence)

    #print(organism_for_tblastn)
    #print(genus)
    #print(taxid)
    #print(title)
    # Display extracted details
    print("\nExtracted Deets:")
    print(f"Protein Title: {selected_info['title']}")
    print(f"Organism: {selected_info['organism']} (Taxonomic ID: {selected_info['taxid']})")
    print(f"Sequence:\n{selected_info['sequence']}")

else:
    print("Failed to retrieve FASTA sequence. Please check the PDB ID.")


def classify_organism_by_taxid(taxid):
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": "taxonomy",
        "id": taxid,
        "retmode": "xml"}
    response = requests.get(url, params=params)
    with open("taxonomy_info.txt","w") as file:
        file.write(response.text)
    if response.status_code != 200:
        print("Error fetching taxonomy data")
        return None
    root = ET.fromstring(response.content)
    #print(root)
    lineage = root.find(".//Lineage").text 
    #print(lineage)
    if "Bacteria" in lineage or "Archaea" in lineage:
        return 9926
    elif "Eukaryota" in lineage:
        return 9927
    else:
        return 1000


classification = classify_organism_by_taxid(taxid)
print(f"Organism with TaxID {taxid} is classified as: {classification}")



#BLASTP
blast_url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
blast_params = {
    "CMD": "Put",
    "PROGRAM": "blastp",
    "DATABASE": "nr",
    "QUERY": sequence,
    "HITLIST_SIZE": 50,  # Request at least 100 results
    "DESCRIPTIONS": 50,  # Ensure 100 descriptions
    "ALIGNMENTS": 50,  # Ensure 100 alignments in results
    "FORMAT_TYPE": "Text",
    "ALIGNMENT_VIEW": "Tabular",  # Set Tabular format
    "OUTFMT":"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp"
}
blast_response = requests.post(blast_url, data=blast_params)
if "RID" in blast_response.text and "RTOE" in blast_response.text:#(Run Time to Estimated Completion)
    rid_start = blast_response.text.find("RID = ") + len("RID = ")
    rid_end = blast_response.text.find("\n", rid_start)
    rid = blast_response.text[rid_start:rid_end].strip()

    rtoe_start = blast_response.text.find("RTOE = ") + len("RTOE = ")
    rtoe_end = blast_response.text.find("\n", rtoe_start)
    rtoe = int(blast_response.text[rtoe_start:rtoe_end].strip())

    print(f"BLAST job submitted successfully. RID: {rid}, Estimated time: {rtoe} seconds.")
else:
    print("Failed to submit BLAST job. Please try again.")
    exit()
    


print("Waiting for search to complete...")
time.sleep(rtoe)
while True:
    result_params = {
        "CMD": "Get",
        "FORMAT_OBJECT": "SearchInfo",
        "RID": rid,
    }
    result_response = requests.get(blast_url, params=result_params)

    if "Status=WAITING" in result_response.text:
        print("Searching...")
        time.sleep(60)
        continue

    if "Status=FAILED" in result_response.text:
        print("Search failed. Please contact NCBI support.")
        exit(4)

    if "Status=UNKNOWN" in result_response.text:
        print("Search expired.")
        exit(3)

    if "Status=READY" in result_response.text:
        if "ThereAreHits=yes" in result_response.text:
            print("Search complete, retrieving results...")
            break
        else:
            print("No hits found.")
            exit(2)
result_params = {
    "CMD": "Get",
    "FORMAT_TYPE": "Text",
    "ALIGNMENT_VIEW": "Tabular",  # Ensure Tabular format
    "RID": rid,
    "OUTFMT": "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp"
}
result_response = requests.get(blast_url, params=result_params)

result_params_text = {
    "CMD": "Get",
    "FORMAT_TYPE": "Text",
    "ALIGNMENT_VIEW": "Text",  # Standard text alignment view
    "RID": rid,
    "OUTFMT": "0"  # Text format
}
result_response_text = requests.get(blast_url, params=result_params_text)

print("\nFull BLAST Results:\n")
print(result_response.text)

output_file = os.path.join(BASE_DIR, "Results", "blast.txt")
with open(output_file, "w",encoding="utf-8") as file:
    file.write(result_response.text)
# Parse and extract accession IDs using regex



import pandas as pd
from collections import defaultdict
import re

def parse_blast_results(file_path):
    """
    Parses a BLAST tabular output file and dynamically determines query length for calculating query coverage.

    :param file_path: Path to the BLAST output file
    :return: Dictionary with accession IDs as keys and query coverage & percentage identity as values
    """

    # Dictionary to store results
    blast_results = defaultdict(lambda: {"query_coverage": 0, "percent_identity": 0, "max_length": 0})#passing lambda to the constructor and using defaultdic cuz i want to dynamically add missing keys on the fly
    
    total_query_length = pdb_length  # To dynamically find query length
    max_qend = 0  # To track the highest query end position if length is not explicitly given

    # Read the file
    with open(file_path, "r") as file:
        in_table = False  # Flag to indicate when we reach the table

        for line in file:
            line = line.strip()

            # Detect start of the table (first actual data row)
            if line.startswith("# Fields:"):
                in_table = True  # Start processing from the next line
                continue  # Restart loop
            
            # Skip all lines until the table starts
            if not in_table or not line: 
                continue
            
            # Ensure the line contains tab-separated values before processing
            columns = line.split("\t")
            if len(columns) < 12:
                continue  
                
            query_id = columns[0]  # Query ID
            accession = columns[1] # Subject Accession ID
            if len(accession)>4 and accession[4]=='_': #we donot want PDB ids since that would give us 100% match
                continue
            percent_identity = float(columns[2])  # % Identity
            alignment_length = int(columns[3])# Alignment Length
            
            query_start = int(columns[6])  # Query Start Position
            query_end = int(columns[7])  # Query End Position
            evalue=float(columns[10])
            # Update max_qend to infer query length if necessary
            max_qend = max(max_qend, query_end)

            # Calculate query coverage (if query length is known)
            if total_query_length is None:
                total_query_length = max_qend  # Default to max qend if query length wasn't found in metadata
            
            query_coverage = ((query_end - query_start + 1) / total_query_length) * 100

            # Update the dictionary
            if accession in blast_results:
                # Accumulate query coverage
                blast_results[accession]["query_coverage"] += query_coverage

                # Keep the highest percent identity from the longest match
                if alignment_length > blast_results[accession]["max_length"]:
                    blast_results[accession]["percent_identity"] = percent_identity
                    blast_results[accession]["max_length"] = alignment_length
            else:
                blast_results[accession] = {
                    "query_coverage": query_coverage,
                    "percent_identity": percent_identity,
                    "max_length": alignment_length,
                    "Evalue":evalue
                }
#Query_4474641	ALB17078.1	100.000	210	0	0	3	212	1	210	3.11e-152	431	100.00
    # Remove "max_length" field from final output
    for accession in blast_results:
        del blast_results[accession]["max_length"]

    return blast_results

# Example Usage
blast_file =os.path.join(BASE_DIR, "Results", "blast.txt") # Replace with the actual file path
blast_data = parse_blast_results(blast_file)

# Convert to DataFrame for easy viewing
df = pd.DataFrame.from_dict(blast_data, orient='index')# from_dict coverts dictionary to pandas dataframe and orient="index" tells pandas to use keys as rows and values(inner dictionaries) as columns
df.to_csv("checking.csv",index=True)
#filtered_accessions = df[(df["query_coverage"] > 85) & (df["percent_identity"] >= 95)].index.tolist()#index means it takes the index row values and tolist means it literally converts to a list
filtered_accessions = df[(df["query_coverage"] > 85) & (df["percent_identity"] >= 95)]
filtered_accessions_list = list(filtered_accessions.itertuples(index=True, name=None))
#print(filtered_accessions_list)
best_accession_from_blastp = min(
    filtered_accessions_list, 
    key=lambda x: (-x[1], -x[2], x[3])  # Sort by max query coverage, max identity, min e-value. x here is each element in the list and min() function does sorting based on minimum hence i do -x[1] and -x[2]
)#also note that this is a custom sorting order
best_acc_bp=best_accession_from_blastp[0]
best_accession_list_bp=list(filtered_accessions_list)

#I'm gonna store the list of best_accessions in a list just in case i want to access them again
final_best_accession_list_bp=[x[0] for x in best_accession_list_bp]
#print(best_acc_bp)
#print(final_best_accession_list_bp)




from Bio import Entrez
from xml.etree import ElementTree


def taxonomy_finder(accession):
    Entrez.email = "your_email@example.com"      
    handle = Entrez.esearch(db="protein", term=accession)
    
    record = Entrez.read(handle)
    handle.close()
    
    uid = record['IdList'][0]
    
    handle = Entrez.efetch(db="protein", id=uid, rettype="gb", retmode="xml")
    records = Entrez.read(handle)
    #print(records)
    handle.close()
    
    
    organism = records[0]['GBSeq_organism']
    features = records[0]['GBSeq_feature-table']
    
    taxid = None
    for feature in features:
        if feature['GBFeature_key'] == 'source':
            for qualifier in feature['GBFeature_quals']:
                if qualifier['GBQualifier_name'] == 'db_xref' and 'taxon:' in qualifier['GBQualifier_value']:
                    taxid = qualifier['GBQualifier_value'].split(':')[1]
                    break
    
    print(f"Organism: {organism}")
    print(f"TaxID: {taxid}")
    return organism,taxid


for accession in final_best_accession_list_bp:
    organism_bp,taxid_bp=taxonomy_finder(accession)
    word=organism_bp.split()
    if len(word)<=1:
        continue
    else:
        break
    
       
counter=1            
if counter==1:
    def get_ncbi_protein_graphics_link(accession):
        base_url = "https://www.ncbi.nlm.nih.gov/protein/"
        link = f"{base_url}{accession}?report=graph"
        return link
 
    link=get_ncbi_protein_graphics_link(best_acc_bp)
    print(f"View the functional annotation of the protein here: {link}")
    

url2 = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
params = {
    "db": "protein",
    "id": best_acc_bp,
    "rettype": "fasta",
    "retmode": "text"
}

response = requests.get(url2, params=params)

if response.status_code == 200:
    fasta_lines = response.text.strip().split("\n")
    print(fasta_lines)
    sequence_bp = ""
    for every_line in fasta_lines:
        if not every_line.startswith(">"):
            sequence_bp+=every_line
                                     
    print(f"Extracted sequence:\n{sequence_bp}")
    
else:
    print(f"Failed to retrieve the sequence. HTTP Status Code: {response.status_code}")

print(organism_bp)


database = "refseq_representative_genomes"
counter_for_tbn=0
RID_for_tbn=""

options = Options()
options.add_experimental_option("prefs", {
    "download.default_directory": "/home/eukprogs/downloads",
    "download.prompt_for_download": False,
    "directory_upgrade": True,
    "safebrowsing.enabled": True
})

# Safe Chrome args for Docker
options.add_argument('--headless=new')  # Headless mode (new version)
options.add_argument('--no-sandbox')  # Required in Docker
options.add_argument('--disable-dev-shm-usage')  # Fixes shared memory crash in Docker

# Use a unique, container-local user data dir
options.add_argument(f"--user-data-dir=/home/eukprogs/chrome_user_data_{os.getpid()}")


driver = webdriver.Chrome(options=options)
#path = r"C:\Users\chira\Desktop\ST_LIFE\query.fasta"
#with open(path, "r") as file:
    #lines = file.readlines()

#protein_seq = "".join(line.strip() for line in lines if not line.startswith(">"))
protein_seq=sequence_bp

def run_tblastn(database,organism,taxid):
    try:
        # Open BLAST page
        driver.get("https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=tblastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome")

        # Wait for the page to load
        wait = WebDriverWait(driver, 10)

        # Click and Enter the Protein Sequence
        seq_box = wait.until(EC.presence_of_element_located((By.ID, "seq")))
        seq_box.click()
        time.sleep(1)
        seq_box.clear()
        seq_box.send_keys(protein_seq)
        time.sleep(1)

        # **Trigger an input event to ensure BLAST recognizes it*
        driver.execute_script("""
            let event = new Event('input', { bubbles: true });
            arguments[0].dispatchEvent(event);
        """, seq_box)
        time.sleep(2)

        
        job_title_box = wait.until(EC.presence_of_element_located((By.ID, "qtitle")))
        job_title_box.clear()
        job_title_box.send_keys("tblastn_eukprogs")

        time.sleep(1)
        # **Ensure Correct Database Selection**
        db_dropdown = wait.until(EC.presence_of_element_located((By.NAME, "DATABASE")))
        db_dropdown.click()
        time.sleep(1)

        # **Use JavaScript to set the correct database option**
        driver.execute_script(f"""
            let dropdown = document.querySelector("[name='DATABASE']");
            dropdown.value = '{database}';
            dropdown.dispatchEvent(new Event('change', {{ bubbles: true }}));
        """)
        print(f"[INFO] Selected {database} database.")

        # **Enter Organism & TaxID, Then Wait for Suggestions**
        organism_box = wait.until(EC.element_to_be_clickable((By.ID, "qorganism")))
        organism_box.clear()  
        organism_box.send_keys(f"{organism} (taxid:{taxid})")
        print(f"[INFO] Entered Organism: {organism} (taxid:{taxid})")
        time.sleep(1)
        # **Check if a matching suggestion appears**
        try:
            matching_suggestion = wait.until(
                EC.presence_of_element_located((By.XPATH, f"//li[contains(text(), '{organism} (taxid:{taxid})')]"))
            )
            matching_suggestion.click()
            print("[INFO] Selected matching organism suggestion.")
        except:
            print("[WARNING] No matching organism suggestion found. Proceeding anyway.")

        # Click the BLAST Button
        time.sleep(2)
        blast_button = wait.until(EC.element_to_be_clickable((By.CLASS_NAME, "blastbutton")))
        blast_button.click()
        print("[INFO] BLAST search submitted successfully!")

        # Wait for BLAST processing
        time.sleep(15)  # i still don't know if i should increase his but this works well so far
        # RID_for_tbn
        global RID_for_tbn
        rid_element = wait.until(EC.presence_of_element_located((By.XPATH, "/html/body/main/section[1]/div[1]/dl/dd[2]/a")))
        RID_for_tbn = rid_element.text
        print(f"[INFO] Extracted RID: {RID_for_tbn}")
        
        
        
        # Check if Alignments Tab is present (indicating results found)
        try:
            alignments_tab = wait.until(
                EC.presence_of_element_located((By.XPATH, "/html/body/main/section[2]/ul/li[3]/button"))
            )
            alignments_tab.click()
            print("[INFO] Opened Alignments tab.")

            # Click the Download button
            download_button = wait.until(
                EC.element_to_be_clickable((By.XPATH, "/html/body/main/section[2]/ul/li[3]/div/div[1]/ul[2]/li/button"))
            )
            download_button.click()
            print("[INFO] Opened Download menu.")

            # Click "Hit Table (text)"
            hit_table_option = wait.until(
                EC.element_to_be_clickable((By.XPATH, "/html/body/main/section[2]/ul/li[3]/div/div[1]/ul[2]/li/ul/li[6]/a"))
            )
            hit_table_option.click()
            print("Hit found, result obtained.")
            return True

        except:
            print("[WARNING] No hits found.")

    except Exception as e:
        print(f"[ERROR] {e}")

    return False

# **First Try with RefSeq**
success = run_tblastn("refseq_representative_genomes", organism_bp, taxid_bp)

# **If No Hits, Retry with Alternatives**
if not success:
    counter_for_tbn
    counter_for_tbn += 1
    success2 = run_tblastn("refseq_representative_genomes", organism_for_tblastn, taxid)
    if not success2:
        success3 = run_tblastn("refseq_genomes", organism_bp, taxid_bp)
        if not success3:
            print("[ERROR] All tBLASTn attempts failed. This organism doesn't have a refseq. Please check the Results folder")
            print("Script stopped executing.")
            exit(1)

# Quit WebDriver
time.sleep(10)
driver.quit()


import shutil
shutil.rmtree(f"/home/eukprogs/chrome_user_data_{os.getpid()}", ignore_errors=True)


tblastn_path = os.path.join(BASE_DIR, "downloads", f"{RID_for_tbn}-Alignment.txt")
gff_file_path = os.path.join(BASE_DIR, "Results", "gff_results.txt")
output_file = "gene_regions.txt"
start_pos_not_found=None
def parse_tblastn_alignment(file_path):
    """
    Extracts the hit sequence ID, start, and end positions from the tBLASTn alignment file.
    """
    with open(file_path, "r") as f:
        lines = f.readlines()
    accession_counter=0
    accession = ""
    start_pos = ""
    end_pos = ""
    start_pos_m3=""
    end_pos_p3=""
    end_pos_m3=""
    for line in lines:
        # Extract Accession ID
        #if "Sequence ID:" in line:
            #accession = line.split(":")
        if "Accession" in line:
            accession_counter+=1
            continue
        if accession_counter==1:
            accession=line.split()[-1]
            accession_counter-=1
            

        # Extract Start & End Positions
        if "Range" in line:
            match = re.search(r"Range\s+\d+:\s+(\d+)\s+to\s+(\d+)", line)
            #if everything fails then please remember to use the next id from the result. 
            if match:
                start_pos = int(match.group(1))
                start_pos_m3=int(match.group(1))-3
                end_pos = int(match.group(2))
                end_pos_p3=int(match.group(2))+3
                end_pos_m3=int(match.group(2))-3
                
                break

    if not accession or not start_pos or not end_pos:
        print("[ERROR] Could not extract accession or positions.")
        return None

    # Ensure correct start and end positions (handle negative strand alignments)
    if start_pos > end_pos:
        start_pos, end_pos = end_pos, start_pos 

    return accession, start_pos, start_pos_m3, end_pos, end_pos_p3, end_pos_m3


def download_gff_file(accession, output_path):
    gff_url = f"https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id={accession}"
    print(f"[INFO] Downloading GFF file for {accession}...")
    
    response = requests.get(gff_url)
    if response.status_code == 200:
        with open(output_path, "w", encoding="utf-8") as f:
            f.write(response.text)
        print(f"[INFO] GFF file saved to {output_path}")
    else:
        print(f"[ERROR] Failed to download GFF file. Status code: {response.status_code}")
        return None

    return output_path



def find_hit(gff_file,start_pos, start_pos_m3, end_pos, end_pos_p3, end_pos_m3):
    line_count=0
    print(start_pos)
    print(end_pos)
    with open(gff_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                line_count+=1
                continue  # Skip comments

            cols = line.strip().split("\t")
            
            if len(cols) < 9:
                line_count+=1
                continue  
            start=int(cols[3])
            end=int(cols[4])
            feature_type=cols[2]
            line_count+=1
            if start == start_pos or end == end_pos_p3:
                if(feature_type=="CDS"):
                    break
            elif start==start_pos_m3 or end==end_pos:
                if(feature_type=="CDS"):
                    break
            elif start==start_pos_m3 or end==end_pos_m3:
                if(feature_type=="CDS"):
                    break
            elif start>start_pos:#start_pos>start esentially means i'm trying to find the closest number to start
                 if(feature_type=="CDS"):
                    line_count+=-1
                    global start_pos_not_found
                    start_pos_not_found=404
                    break
    return line_count

def find_upstr_dowstr(gff_file,line_count):
    upstream_count = max(1, line_count - 5) 
    downstream_count = line_count + 5
    counter = 0
    genes_dict = {}

    with open(gff_file, "r") as f:
        for line in f:
            counter += 1
            if counter < upstream_count or counter > downstream_count:
                continue  #essentially continue will skip any remaining code and moves the cursor to the next line in the next iteration

            cols = line.strip().split("\t")
            if len(cols) < 9:
                continue
            
            
            #taking only product and nothing else for now from the attributes field
            if start_pos_not_found==404:
                print("There was no matching start position, taking the closest start position from the GFF file")
                attributes=cols[8]
                match = re.search(r'product=([^;]+)', attributes)
                if cols[2] == "CDS" or cols[2] == "exon":
                    genes_dict[counter] = {
                    "Start": cols[3],
                    "End": cols[4],
                    "Strand": cols[6],
                    "Attributes": match.group(1)
                    }
            else:          
                attributes=cols[8]
                match = re.search(r'product=([^;]+)', attributes)
                if cols[2] == "CDS" or cols[2] == "exon":
                    genes_dict[counter] = {
                        "Start": cols[3],
                        "End": cols[4],
                        "Strand": cols[6],
                        "Attributes": match.group(1)
                    }

    return genes_dict   
    
alignment_info = parse_tblastn_alignment(tblastn_path)
if alignment_info:
    accession, start_pos, start_pos_m3, end_pos, end_pos_p3, end_pos_m3 = alignment_info

final_refseq_accession=alignment_info[0]

output_path=download_gff_file(accession,gff_file_path)#a variable to hold the path since the function has to return something, but i won't use this anywhere


if gff_file_path:
    line_count = find_hit(gff_file_path,start_pos, start_pos_m3, end_pos, end_pos_p3, end_pos_m3)
    stream=find_upstr_dowstr(gff_file_path,line_count)




labels = ["A", "B", "C", "D", "E"]
lines = []
start_values = []
counter=0
"""
HOW TO USE ENUMERATE:
>>> list1=["a","b","c","d"]
>>> print(enumerate(list1))
<enumerate object at 0x000002171FF94B80>
>>> print(list(enumerate(list1)))
[(0, 'a'), (1, 'b'), (2, 'c'), (3, 'd')]


HOW IDX THING WORKS:
idx → The index of the current loop iteration (0, 1, 2...).
key → The dictionary key (e.g., 998, 1000).
value → The dictionary value (which is another dictionary in this case) which is our original target. the keys are just line numbers which are useless here
"""
for idx, (key, value) in enumerate(stream.items()):# the item()function will return a view object which returns a list of tuples. the tuples contain key-value pairs of a dictionary
    attribute = value['Attributes']
    counter+=1
    if counter==3:
        genome_annot=attribute
    strand = value['Strand']

    start_values.append(value['Start'])  #OR  just use start=1 in the enumerate function
    lines.append(f"{labels[idx]}: {attribute}({strand})") #since labels is a list and follows zero-based indexing, labels[idx] means access zeroth position of labels


start_values_str = f"start_values={','.join(start_values)}"
title_str= f"PDB: {title}"
Genome_accession=f"Genome_Accession: {final_refseq_accession}"
organism=f"Organism: {organism_for_tblastn}"
# Write formatted content to a file
file_path = os.path.join(BASE_DIR, "Results", "stream.txt") # Change this path as needed
with open(file_path, "w") as f:
    f.write("\n".join(lines) + "\n" + "\n") #\n".join(lines)-means joining a list of strings with a new line character
    f.write(start_values_str + "\n" + "\n")
    f.write(Genome_accession+ "\n" + "\n")
    f.write(organism+ "\n" + "\n")
    f.write(title_str+ "\n" + "\n")
    

print(f"[INFO] Data successfully written") 
import subprocess
#from IPython.display import Image, display

# Path to your R script
r_script_path = os.path.join(BASE_DIR,"ggplot.R")

# Run the R script
subprocess.run(["Rscript", r_script_path], check=True)

# Display the generated PNG (assumes the script outputs this exact file)
output_image_path = os.path.join(BASE_DIR,"Results","genomic_annotation_plot.png")

# Display the image inside Jupyter Notebook
#display(Image(output_image_path))



def gemini_api(api_key,title,gennome_annot,selected_chain_group,final_refseq_accession):

	model_list = ["models/gemini-1.5-pro-latest"]
	model = genai.GenerativeModel(model_list[0])
	
	genai.configure(api_key=api_key)
	response = model.generate_content(
		     
		    f"""Determine if the following two protein annotations describe the same biological function.

		    PDB: {title}  
		    Genome: {genome_annot}
		     	
		    Instructions:
		    1. Ignore keyword similarity; focus on actual biological roles.  
		    2. Consider two annotations "consistent" if:
	   	    - They describe the same enzyme class/function.
	   	    - Both are uncharacterized/hypothetical/DUF proteins.
	   	    - The difference is wording but the role is the same (e.g. tRNA-intron lyase vs tRNA-splicing endonuclease).
		    3. Otherwise, label "inconsistent" if roles differ clearly (e.g. synthetase vs transcription factor).

		    Answer strictly "consistent" or "inconsistent":""",
		    generation_config={"temperature": 0.0}
		)
	answer = response.text.strip()
	print(f"PDB title:{title}")
	print(f"Gene annotation:{genome_annot}")
	print(f"Flag:{response.text}")
	
	output_csv =os.path.join(BASE_DIR,"Results","Final_results.csv")
	with open(output_csv, mode="w", newline="") as file:
    		writer = csv.writer(file)
    		writer.writerow(["pdb_id", "genome_accession", "pdb_description", "gene_annotation", "flag"])
    		writer.writerow([selected_chain_group, final_refseq_accession, title, genome_annot, response.text])
	
	print("Generative AI is experimental, please do not make any decision based on the classification the AI does")	 






if api_key:
	gemini_api(api_key,title,gennome_annot,selected_chain_group,final_refseq_accession)

else:
	with open(output_csv, mode="w", newline="") as file:
    		writer = csv.writer(file)
    		writer.writerow(["pdb_id", "genome_accession", "pdb_description", "gene_annotation"])
    		writer.writerow([selected_chain_group, final_refseq_accession, title, genome_annot])
    			















