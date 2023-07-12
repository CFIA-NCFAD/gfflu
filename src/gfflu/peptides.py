from pathlib import Path

iav_faa = Path(__file__).parent.joinpath("data", "iav-annotation.faa")


class Peptides:
    """Short peptide reference sequences for Influenza A virus.

    Attributes
    ----------
    ha_signal_peptides : tuple
        Peptide sequences for the HA signal peptide region of gene segment 4 from the
        `NCBI Influenza Virus Resource segment 4 sequences
        <https://ftp.ncbi.nih.gov/genomes/INFLUENZA/ANNOTATION/PROTEIN-A/Seg4/>`_.
        The original sequence ID (e.g. "seg4mature3A1") is retained.

    pb1_f2_peptides : tuple
        Peptide sequences for the PB2-F2 gene of IAV gene segment 2 from the
        `NCBI Influenza Virus Resource segment 2 sequences
        <https://ftp.ncbi.nih.gov/genomes/INFLUENZA/ANNOTATION/PROTEIN-A/Seg2/>`_.
        The original sequence ID (e.g. "seg2prot2A") is retained.
    """

    ha_signal_peptides = (
        ("HA|signal_peptide_region_of_CDS||seg4mature3A1", "MKTIIVLSYLFCLALS"),
        ("HA|signal_peptide_region_of_CDS||seg4mature5A1", "MEKIVLLLAIVSLVKS"),
        ("HA|signal_peptide_region_of_CDS||seg4mature5B1", "MEKIVLLLAIVSLVKS"),
        ("HA|signal_peptide_region_of_CDS||seg4mature7A1", "MNIQILAFIACVLTGAKG"),
        ("HA|signal_peptide_region_of_CDS||seg4mature7C1", "MNIQILVAIACALIETKA"),
        ("HA|signal_peptide_region_of_CDS||seg4mature7D1", "MNTQILVFALVAIIPTNA"),
        ("HA|signal_peptide_region_of_CDS||seg4mature7E1", "MNTQILILTLVAAIHTNA"),
        ("HA|signal_peptide_region_of_CDS||seg4mature7F1", "MNTQILVFALVAIIPTNA"),
        ("HA|signal_peptide_region_of_CDS||seg4matureA1", "MKTILVVLLYTFTTTYA"),
        ("HA|signal_peptide_region_of_CDS||seg4matureB1", "MKTIIAFSYILCLIFA"),
        ("HA|signal_peptide_region_of_CDS||seg4matureC1", "MKTTIILILLTHWVYS"),
        ("HA|signal_peptide_region_of_CDS||seg4matureD1", "MAIIYLILLFTAVRG"),
        ("HA|signal_peptide_region_of_CDS||seg4matureE1", "MIALILVALALSHTAYS"),
        ("HA|signal_peptide_region_of_CDS||seg4matureF1", "MIAIIVVAILATAGRS"),
        ("HA|signal_peptide_region_of_CDS||seg4matureG1", "MLSIVILFLLIAENSS"),
        ("HA|signal_peptide_region_of_CDS||seg4matureH1", "MALNVIATLTLISVCVHA"),
        ("HA|signal_peptide_region_of_CDS||seg4matureI1", "MEKFIAIATLASTNAY"),
        ("HA|signal_peptide_region_of_CDS||seg4matureJ1", "METKAIIAALLMVTAANA"),
        ("HA|signal_peptide_region_of_CDS||seg4matureK1", "MEKTLLFAAIFLCVKA"),
        ("HA|signal_peptide_region_of_CDS||seg4matureL1", "MEKFIILSTVLAASFAY"),
        ("HA|signal_peptide_region_of_CDS||seg4matureM1", "MVIKVLYFLIVLLSRYSKA"),
        ("HA|signal_peptide_region_of_CDS||seg4matureN1", "MYKVVVIIALLGAVRG"),
        ("HA|signal_peptide_region_of_CDS||seg4matureO1", "MNTQILVFALVAIIPTNA"),
        ("HA|signal_peptide_region_of_CDS||seg4matureP1", "MEKIVLLLAIVSLVKS"),
    )
    pb1_f2_peptides = (
        (
            "PB1-F2|CDS|PB1-F2_protein|seg2prot2A",
            "MEQEQGTPWTQSTEHTNIQKRGSGRQIQKLGHPNSTQLMDHYLRIMSQVDMHKQTVSWRLWPSLKNPTQGSLRTHALKQWKSFNKQGWTN",
        ),
        (
            "PB1-F2|CDS|PB1-F2_protein|seg2prot2B",
            "MEQEQGTPWTQSTEHTNIQKRGSGRQIQKLGHPNSTQLMDHYLRIMSQVDMHKQTVSWRLWPSLKNPTQGSLRTHALKQWKSFNKQG",
        ),
        (
            "PB1-F2|CDS|PB1-F2_protein|seg2prot2C",
            "MEQEQGTPWTQSTEHTNIQKRGSGRQIQKLGHPNSTQLMDHYLRIMSQVDMHKQTVSWRLWPSLKNPTQGSLRTHALKQ",
        ),
        ("PB1-F2|CDS|PB1-F2_protein|seg2prot2D", "MEQEQGTPWTQSTEHTNIQKRGSGRQIQKLGHPNSTQLMDHYLRIMSQVDMHKQTVS"),
        ("PB1-F2|CDS|PB1-F2_protein|seg2prot2E", "MEQEQGTPWTQSTEHTNIQKRGSGRQIQKLGHPNSTQLMDHYLRIMSQVDMHKQTVSWR"),
        (
            "PB1-F2|CDS|PB1-F2_protein|seg2prot2F",
            "MEQEQGTPWTQSTEHTNIQKRGSGRQIQKLGHPNSTQLMDHYLRIMSQVDMHKQTVSWRLWPSLKNPTQGSLRTHALKQWKSFNKQGWTNLLKVARLMIGH",
        ),
        ("PB1-F2|CDS|PB1-F2_protein|seg2prot2G", "MEQEQGTPWTQSTEHTNIQKRGSGRQIQKLGHPNSTQLMDHYLRIMSQVDMHKQTVSWRLWPS"),
        ("PB1-F2|CDS|PB1-F2_protein|seg2prot2H_putative", "MDHYLRIMSQVDMHKQTVSWRLWPSLKNPTQGSLRTHALKQWKSFNKQGWTN"),
        (
            "PB1-F2|CDS|PB1-F2_protein|seg2prot2I",
            "MEQEQGTPWTQSTEHTNIQKRGSGRQIQKLGHPNSTQLMDHYLRIMSQVDMHKQTVSWRLWPSLKNPTQGSLRTHALK",
        ),
        ("PB1-F2|CDS|PB1-F2_protein|seg2prot2J", "MEQEQGTPWTQSTEHTNIQKRGSGRQIQKLGHPNSTQLMDHYLRIMSQVDMHKQTV"),
        (
            "PB1-F2|CDS|PB1-F2_protein|seg2prot2K",
            "MEQEQGTPWTQSTEHTNIQKRGSGRQIQKLGHPNSTQLMDHYLRIMSQVDMHKQTVSWRLWPSLKNPTQGSLRTHALKQWR",
        ),
        (
            "PB1-F2|CDS|PB1-F2_protein|seg2prot2L",
            "MEQEQGTPWTQSTEHTNIQKRGSGRQIQKLGHPNSTQLMDHYLRIMSQVDMHKQTVSWRLWPSLKNPTQGSLRTHA",
        ),
        (
            "PB1-F2|CDS|PB1-F2_protein|seg2prot2M",
            "MEQEQDTPWTQSTEHTNIQKKGNGRQIQRLGHPSSIRLMDHYLKIMNQVAMHKQTVSWRPWLSLKNPTQGYLRIHALKQWKLSNKQGWIN",
        ),
    )
