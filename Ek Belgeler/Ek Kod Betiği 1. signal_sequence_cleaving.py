
#CODE CLEAVES THE POSITION AND THE REMOVED RESIDUES NUMBER IS SAME WITH THE CLEAVAGE SITE POSITION

#!/usr/bin/env python3
"""
Protein Cleavage Site Analysis Script
Processes protein sequences by finding cleavage sites and removing the N-terminal region.
Outputs cleaved sequences to a FASTA file.
"""

def find_cleavage_and_cut(sequence, cleavage_position):
    """
    Cut the sequence at the specified position, removing everything before it.
    
    Args:
        sequence (str): Protein sequence
        cleavage_position (int or str): Position in the sequence to cleave (1-based indexing)
    
    Returns:
        str: Cleaved sequence (from cleavage position onwards)
    """
    # Convert to integer if needed
    try:
        cleavage_position = int(cleavage_position)
    except (ValueError, TypeError):
        print(f"Error: Cleavage position '{cleavage_position}' is not a valid number")
        return sequence
    
    # Convert to 0-based indexing (cleavage happens after the specified position)
    position = cleavage_position
    
    if position < 0 or position >= len(sequence):
        print(f"Warning: Cleavage position {cleavage_position} is out of range for sequence length {len(sequence)}")
        return sequence  # Return original if position is invalid
    
    # Return sequence starting from the cleavage position
    cleaved_seq = sequence[position:]
    return cleaved_seq


def write_fasta(output_file, sequences_dict):
    """
    Write sequences to a FASTA file.
    
    Args:
        output_file (str): Path to output FASTA file
        sequences_dict (dict): Dictionary with protein names as keys and sequences as values
    """
    with open(output_file, 'w') as f:
        for protein_name, sequence in sequences_dict.items():
            # Write FASTA header
            f.write(f">{protein_name}\n")
            
            # Write sequence in 60-character lines (standard FASTA format)
            for i in range(0, len(sequence), 60):
                f.write(sequence[i:i+60] + "\n")


def process_proteins(protein_sequences, cleavage_sites, output_file="cleaved_proteins.fasta"):
    """
    Main function to process proteins with their cleavage sites.
    
    Args:
        protein_sequences (dict): Dictionary of protein names to sequences
        cleavage_sites (dict): Dictionary of protein names to cleavage positions (integers)
        output_file (str): Output FASTA file name
    """
    cleaved_sequences = {}
    
    # Process each protein
    for protein_name in protein_sequences:
        if protein_name in cleavage_sites:
            sequence = protein_sequences[protein_name]
            # Remove any whitespace from the sequence
            sequence = ''.join(sequence.split())
            
            cleavage_position = cleavage_sites[protein_name]
            
            print(f"Processing {protein_name}...")
            print(f"  Original length: {len(sequence)}")
            print(f"  Cleavage position: {cleavage_position}")
            
            # Find cleavage site and cut
            cleaved_seq = find_cleavage_and_cut(sequence, cleavage_position)
            
            print(f"  Cleaved length: {len(cleaved_seq)}")
            print(f"  Residues removed: {len(sequence) - len(cleaved_seq)}")
            
            cleaved_sequences[protein_name] = cleaved_seq
        else:
            print(f"Warning: No cleavage site found for {protein_name}")
    
    # Write to FASTA file
    if cleaved_sequences:
        write_fasta(output_file, cleaved_sequences)
        print(f"\nCleaved sequences written to {output_file}")
        print(f"Total proteins processed: {len(cleaved_sequences)}")
    else:
        print("No sequences to write!")

# Example usage
if __name__ == "__main__":
    # Example protein sequences dictionary
    protein_sequences = {
        "BaCut1": "MKANLILACASVVAGVSAAPLERRDGACSKYTIIDTRGTGELQGPSTGFLTMNRNILSEVPGGDQYHTVY PADWSQISTQGTLDIVNKVYSTLESDPNHCFVLEGYSQGAAATVNALPKLTGDSFDAVKAVFLIGNPMHK SGLECNVDTLGGKSTAYANGLSVVLGSIPDEWVSKTMDVCNFGDGVCDTLTGVGITLQHLDYPLDANVQN MGTEFVVKALTS",
        "CpCut1": "MKFTTAVALLAGTAVALPTAVERRQSVGSTVNELVSGNCRPITFIFARGSTEIGNIGSSVGPPTCQGLKS TYGANNVACQGVGSPYDATLAANLLPEGTTSAAYGEAQRLFNLANTKCPNTKIVAGGYSQGAAVMVASVR RLSTAVKNKIAGVVLYGNTRNAQEGGKIPSFPPEKALTICNASDGVCGGGLVVTPGHLTYQVDVPKAVDY LESRIAAAGGI",
        "HiC": "MQLGAIENGLESGSANACPDAILIFARGSTEPGNMGITVGPALANGLESHIRNIWIQGVGGPYDAALATN FLPRGTSQANIDEGKRLFALANQKCPNTPVVAGGYSQGAALIAAAVSELSGAVKEQVKGVALFGYTQNLQ NRGGIPNYPRERTKVFCNVGDAVCTGTLIITPAHLSYTIEARGEAARFLRDRIRA",
        "LCC": "MDGVLWRVRTAALMAALLALAAWALVWASPSVEAQSNPYQRGPNPTRSALTADGPFSVATYTVSRLSVSG FGGGVIYYPTGTSLTFGGIAMSPGYTADASSLAWLGRRLASHGFVVLVINTNSRFDYPDSRASQLSAALN YLRTSSPSAVRARLDANRLAVAGHSMGGGGTLRIAEQNPSLKAAVPLTPWHTDKTFNTSVPVLIVGAEAD TVAPVSQHAIPFYQNLPSTTPKVYVELDNASHFAPNSNNAAISVYTISWMKLWVDNDTRYRQFLCNVNDP ALSDFRTNNRHCQ",
        "PudA": "MNSRSLTKAIRFPTILALAGFSVLGACGGSDNDSSSNNQGAPAVAITVAGQVQAVDRLGMRRYFGIPFAA PPVGNLRWMPPAPPQSWAAPLAKTQSNAPCMQTGATDPLRLPNGTEDCLYLDVHAPATGEGPFPVMVWIH GGAFSIGGTITYADPSPLVSKGVIVVNIAYRMGAMGFLGHPSLRAADGTVGNYGIMDQQAALRWVQDNIA AFGGDKSNVTIFGESAGGFSVMTHLASPLSKGLFAKAIVQSGGYGFDRQLTQAQLEAQSTSIVNSALAAA GVSCPTVDAACLRGLSAELVNNQLATAFTTANWSPVPSVDGKVLPKSIKATFVAGENNKVPLVNGSNQDE WSYFVASRELVAGPLTAAQYPSYLQTSLGLPPSLATVYPLTDYGTNTAQQPSLAATAAGTDMHFSCPALN LSKRVLSQATPIFMYEFRDRTAIPSIGRNTISFNQGAGHTYELQYLFNLRDLETAEHRDLQASMARYWTN FARTSNPNNGDPVATSWPAFTGPTKVLGLDVASAGGIRELATFETDHKCNTAWTSLTF",
        "PurH": "MKLSRLIAMALAGAALTAGTAIQPVVAAENPYERGPAPTNSSIEATRGPYAVSTKTISSLSARGFGGGTIYYPTSTADGTFGVVAISPGYTAAQSTIQWLGPRIASQGFVVITIDTNTRLDQPGSRGTQLLAALDQTIADTTVRSRIDASRQAVVGHSMGGGGTLEAAKSRRSIEATVGLTPWNLDKTWPEVEAASLEIGAQNDTVAPPGSHAIPFYNSLTNAERRAYLELRGASHFAPNTSNTTIAKYTIAWLKRYVDDDTRYEQFISPGPSPSLTNGISDYRIQ",
        "Tcur0390": "MKRTLKRALSLLPAAALAASALVAASPAQAAANPYQRGPNPTEASITAARGPFNTAEITVSRLSVSGFGGGKIYYPTTTSEGTFGAIAISPGFTAYWSSLEWLGHRLASQGFVVIGIETNTTLDQPDQRGQQLLAALDYLTQRSAVRDRVDASRLAVAGHSMGGGGSLEAAKARTSLKAAIPLAPWNLDKTWPEVRTPTLIIGGELDAVAPVATHSIPFYNSLSNAPEKAYLELDNASHFFPNITNTQMAKYMIAWMKRFIDDDTRYTQFLCPPPSTGLLSDFSDARFTCPM",
        "Tcur1278": "MSLRKSFGLLSATAALVAGLVAAPPAQAAANPYQRGPDPTESLLRAARGPFAVSEQSVSRLSVSGFGGGRIYYPTTTSQGTFGAIAISPGFTASWSSLAWLGPRLASHGFVVIGIETNTRLDQPDSRGRQLLAALDYLTQRSSVRNRVDASRLAVAGHSMGGGGTLEAAKSRTSLKAAIPIAPWNLDKTWPEVRTPTLIIGGELDSIAPVATHSIPFYNSLTNAREKAYLELNNASHFFPQFSNDTMAKFMISWMKRFIDDDTRYDQFLCPPPRAIGDISDYRDTCPHT",
        "TfCut2": "MAVMTPRRERSSLLSRALQVTAAAATALVTAVSLAAPAHAANPYERGPNPTDALLEASSGPFSVSEENVSRLSASGFGGGTIYYPRENNTYGAVAISPGYTGTEASIAWLGERIASHGFVVITIDTITTLDQPDSRAEQLNAALNHMINRASSTVRSRIDSSRLAVMGHSMGGGGTLRLASQRPDLKAAIPLTPWHLNKNWSSVTVPTLIIGADLDTIAPVATHAKPFYNSLPSSISKAYLELDGATHFAPNIPNKIIGKYSVAWLKRFVDNDTRYTQFLCPGPRDGLFGEVEEYRSTCPF",
        "CALA": "MRVSLRSITSLLAAATAAVLAAPATETLDRRAALPNPYDDPFYTTPSNIGTFAKGQVIQSRKVPTDIGNANNAASFQLQYRTTNTQNEAVADVATVWIPAKPASPPKIFSYQVYEDATALDCAPSYSYLTGLDQPNKVTAVLDTPIIIGWALQQGYYVVSSDHEGFKAAFIAGYEEGMAILDGIRALKNYQNLPSDSKVALEGYSGGAHATVWATSLADSYAPELNIVGASHGGTPVSAKDTFTFLNGGPFAGFALAGVSGLSLAHPDMESFIEARLNAKGQQTLKQIRGRGFCLPQVVLTYPFLNVFSLVNDTNLLNEAPIAGILKQETVVQAEASYTVSVPKFPRFIWHAIPDEIVPYQPAATYVKEQCAKGANINFSPYPIAEHLTAEIFGLVPSLWFIKQAFDGTTPKVICGTPIPAIAGITTPSADQVLGSDLANQLRSLNGKQSAFGKPFGPITPP"
    }
    
    # Example cleavage sites dictionary (numeric positions, 1-based)
    cleavage_sites = {
        "BaCut1": 24,
        "CpCut1": 24,
        "HiC": 17,
        "LCC": 34,
        "PudA": 26,
        "PurH": 27,
        "Tcur0390": 30,
        "Tcur1278": 28,
        "TfCut2": 40,
        "CALA": 31
            }
    
    # Process proteins and generate FASTA file
    process_proteins(protein_sequences, cleavage_sites, "cleaved_proteinsall.fasta")