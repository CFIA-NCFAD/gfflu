from pathlib import Path
from tempfile import TemporaryDirectory

from gfflu.annotation import run_annotation

datadir = Path(__file__).parent / "data"


def test_annotate_seg1():
    input_fasta = datadir / "Segment_1_PB2.MH201221.fasta"
    with TemporaryDirectory(prefix="gfflu-") as tmpdir:
        gffrec = run_annotation(input_fasta, Path(tmpdir), prefix=input_fasta.stem)
        assert len(gffrec.features) == 1
        feature = gffrec.features[0]
        assert feature.type == "gene"
        assert feature.qualifiers["gene"] == ["PB2"]
        # zero-based indexing for start using ExactPosition
        assert feature.location.start == 15
        assert feature.location.end == 2295
        assert len(feature.sub_features) == 1
        subfeature = feature.sub_features[0]
        assert subfeature.type == "CDS"
        assert subfeature.qualifiers["gene"] == ["PB2"]
        assert subfeature.qualifiers["product"] == ["polymerase PB2"]
        # check that subfeature location is correct
        assert subfeature.location.start == 15
        assert subfeature.location.end == 2295


def test_annotate_seg2():
    input_fasta = datadir / "Segment_2_PB1.CY147460.fasta"
    with TemporaryDirectory(prefix="gfflu-") as tmpdir:
        gffrec = run_annotation(input_fasta, Path(tmpdir), prefix=input_fasta.stem)
        assert len(gffrec.features) == 2
        feature_pb1 = gffrec.features[0]
        feature_pb1_f2 = gffrec.features[1]
        assert feature_pb1.type == "gene"
        assert feature_pb1.qualifiers["gene"] == ["PB1"]
        # zero-based indexing for start using ExactPosition
        assert feature_pb1.location.start == 12
        assert feature_pb1.location.end == 2286
        assert len(feature_pb1.sub_features) == 1
        subfeature_pb1 = feature_pb1.sub_features[0]
        assert subfeature_pb1.type == "CDS"
        assert subfeature_pb1.qualifiers["gene"] == ["PB1"]
        assert subfeature_pb1.qualifiers["product"] == ["polymerase PB1"]
        # check that subfeature location is correct
        assert subfeature_pb1.location.start == 12
        assert subfeature_pb1.location.end == 2286
        assert feature_pb1_f2.type == "gene"
        assert feature_pb1_f2.qualifiers["gene"] == ["PB1-F2"]
        # zero-based indexing for start using ExactPosition
        assert feature_pb1_f2.location.start == 106
        assert feature_pb1_f2.location.end == 370
        assert len(feature_pb1_f2.sub_features) == 1
        subfeature_pb1_f2 = feature_pb1_f2.sub_features[0]
        assert subfeature_pb1_f2.type == "CDS"
        assert subfeature_pb1_f2.qualifiers["gene"] == ["PB1-F2"]
        assert subfeature_pb1_f2.qualifiers["product"] == ["PB1-F2 protein"]
        # check that subfeature location is correct
        assert subfeature_pb1_f2.location.start == 106
        assert subfeature_pb1_f2.location.end == 370


def test_annotate_seg3():
    input_fasta = datadir / "Segment_3_PA.CY146806.fasta"
    with TemporaryDirectory(prefix="gfflu-") as tmpdir:
        gffrec = run_annotation(input_fasta, Path(tmpdir), prefix=input_fasta.stem)
        assert len(gffrec.features) == 2
        feature_pa = gffrec.features[0]
        feature_pa_x = gffrec.features[1]
        assert feature_pa.type == "gene"
        assert feature_pa.qualifiers["gene"] == ["PA"]
        # zero-based indexing for start using ExactPosition
        assert feature_pa.location.start == 12
        assert feature_pa.location.end == 2163
        assert len(feature_pa.sub_features) == 1
        subfeature_pa = feature_pa.sub_features[0]
        assert subfeature_pa.type == "CDS"
        assert subfeature_pa.qualifiers["gene"] == ["PA"]
        assert subfeature_pa.qualifiers["product"] == ["polymerase PA"]
        # check that subfeature location is correct
        assert subfeature_pa.location.start == 12
        assert subfeature_pa.location.end == 2163
        assert feature_pa_x.type == "gene"
        assert feature_pa_x.qualifiers["gene"] == ["PA-X"]
        # zero-based indexing for start using ExactPosition
        assert feature_pa_x.location.start == 12
        assert feature_pa_x.location.end == 772
        assert len(feature_pa_x.sub_features) == 1
        subfeature_pa_x = feature_pa_x.sub_features[0]
        assert subfeature_pa_x.type == "CDS"
        assert subfeature_pa_x.qualifiers["gene"] == ["PA-X"]
        assert subfeature_pa_x.qualifiers["product"] == ["PA-X protein"]
        # check that subfeature location is correct
        assert subfeature_pa_x.location.start == 12
        assert subfeature_pa_x.location.end == 772
        assert subfeature_pa_x.qualifiers["Frameshift"] == ["1"]


def test_annotation_seg4():
    input_fasta = datadir / "Segment_4_HA.MH201222.fasta"
    with TemporaryDirectory(prefix="gfflu-") as tmpdir:
        gffrec = run_annotation(input_fasta, Path(tmpdir), prefix=input_fasta.stem)
        assert len(gffrec.features) == 1
        feature_ha = gffrec.features[0]
        assert feature_ha.type == "gene"
        assert feature_ha.qualifiers["gene"] == ["HA"]
        # zero-based indexing for start using ExactPosition
        assert feature_ha.location.start == 20
        assert feature_ha.location.end == 1721
        assert len(feature_ha.sub_features) == 4
        subfeature_ha = feature_ha.sub_features[0]
        subfeature_sig_peptide = feature_ha.sub_features[1]
        subfeature_ha1 = feature_ha.sub_features[2]
        subfeature_ha2 = feature_ha.sub_features[3]
        assert subfeature_ha.type == "CDS"
        assert subfeature_ha.qualifiers["gene"] == ["HA"]
        assert subfeature_ha.qualifiers["product"] == ["hemagglutinin"]
        # check that subfeature location is correct
        assert subfeature_ha.location.start == 20
        assert subfeature_ha.location.end == 1721

        assert subfeature_sig_peptide.type == "signal_peptide_region_of_CDS"
        assert subfeature_sig_peptide.location.start == 20
        assert subfeature_sig_peptide.location.end == 71

        assert subfeature_ha1.qualifiers["product"] == ["HA1"]
        assert subfeature_ha1.type == "mature_protein_region_of_CDS"
        assert subfeature_ha1.location.start == 71
        assert subfeature_ha1.location.end == 1052

        assert subfeature_ha2.qualifiers["product"] == ["HA2"]
        assert subfeature_ha2.type == "mature_protein_region_of_CDS"
        assert subfeature_ha2.location.start == 1052
        assert subfeature_ha2.location.end == 1718


def test_annotation_seg5():
    input_fasta = datadir / "Segment_5_NP.MH085254.fasta"
    with TemporaryDirectory(prefix="gfflu-") as tmpdir:
        gffrec = run_annotation(input_fasta, Path(tmpdir), prefix=input_fasta.stem)
        assert len(gffrec.features) == 1
        feature_np = gffrec.features[0]
        assert feature_np.type == "gene"
        assert feature_np.qualifiers["gene"] == ["NP"]
        # zero-based indexing for start using ExactPosition
        assert feature_np.location.start == 43
        assert feature_np.location.end == 1540
        assert len(feature_np.sub_features) == 1
        subfeature_np = feature_np.sub_features[0]
        assert subfeature_np.type == "CDS"
        assert subfeature_np.qualifiers["gene"] == ["NP"]
        assert subfeature_np.qualifiers["product"] == ["nucleocapsid protein"]
        # check that subfeature location is correct
        assert subfeature_np.location.start == 43
        assert subfeature_np.location.end == 1540


def test_annotation_seg6():
    input_fasta = datadir / "Segment_6_NA.EF190976.fasta"
    with TemporaryDirectory(prefix="gfflu-") as tmpdir:
        gffrec = run_annotation(input_fasta, Path(tmpdir), prefix=input_fasta.stem)
        assert len(gffrec.features) == 1
        feature_na = gffrec.features[0]
        assert feature_na.type == "gene"
        assert feature_na.qualifiers["gene"] == ["NA"]
        # zero-based indexing for start using ExactPosition
        assert feature_na.location.start == 20
        assert feature_na.location.end == 1385
        assert len(feature_na.sub_features) == 1
        subfeature_na = feature_na.sub_features[0]
        assert subfeature_na.type == "CDS"
        assert subfeature_na.qualifiers["gene"] == ["NA"]
        assert subfeature_na.qualifiers["product"] == ["neuraminidase"]
        # check that subfeature location is correct
        assert subfeature_na.location.start == 20
        assert subfeature_na.location.end == 1385


def test_annotation_seg7():
    input_fasta = datadir / "Segment_7_M.MH085255.fasta"
    with TemporaryDirectory(prefix="gfflu-") as tmpdir:
        gffrec = run_annotation(input_fasta, Path(tmpdir), prefix=input_fasta.stem)
        assert len(gffrec.features) == 2
        feature_m1 = gffrec.features[0]
        assert feature_m1.type == "gene"
        assert feature_m1.qualifiers["gene"] == ["M1"]
        # zero-based indexing for start using ExactPosition
        assert feature_m1.location.start == 23
        assert feature_m1.location.end == 782
        assert len(feature_m1.sub_features) == 1
        subfeature_m1 = feature_m1.sub_features[0]
        assert subfeature_m1.type == "CDS"
        assert subfeature_m1.qualifiers["gene"] == ["M1"]
        assert subfeature_m1.qualifiers["product"] == ["matrix protein 1"]
        # check that subfeature location is correct
        assert subfeature_m1.location.start == 23
        assert subfeature_m1.location.end == 782

        feature_m2 = gffrec.features[1]
        assert feature_m2.type == "gene"
        assert feature_m2.qualifiers["gene"] == ["M2"]
        # zero-based indexing for start using ExactPosition
        assert feature_m2.location.start == 23
        assert feature_m2.location.end == 1005
        assert len(feature_m2.sub_features) == 2
        subfeature_m2_1 = feature_m2.sub_features[0]
        subfeature_m2_2 = feature_m2.sub_features[1]
        assert subfeature_m2_1.type == "CDS"
        assert subfeature_m2_1.qualifiers["gene"] == ["M2"]
        assert subfeature_m2_1.qualifiers["product"] == ["matrix protein 2"]
        # check that subfeature location is correct
        assert subfeature_m2_1.location.start == 23
        assert subfeature_m2_1.location.end == 49

        assert subfeature_m2_2.type == "CDS"
        assert subfeature_m2_2.qualifiers["gene"] == ["M2"]
        assert subfeature_m2_2.qualifiers["product"] == ["matrix protein 2"]
        # check that subfeature location is correct
        assert subfeature_m2_2.location.start == 737
        assert subfeature_m2_2.location.end == 1005


def test_annotation_seg8():
    input_fasta = datadir / "Segment_8_NS.MH085256.fasta"
    with TemporaryDirectory(prefix="gfflu-") as tmpdir:
        gffrec = run_annotation(input_fasta, Path(tmpdir), prefix=input_fasta.stem)
        assert len(gffrec.features) == 2
        if gffrec.features[0].qualifiers["gene"] == ["NS1"]:
            feature_ns1 = gffrec.features[0]
            feature_ns2 = gffrec.features[1]
        assert feature_ns1.type == "gene"
        # zero-based indexing for start using ExactPosition
        assert feature_ns1.location.start == 24
        assert feature_ns1.location.end == 717
        assert len(feature_ns1.sub_features) == 1
        subfeature_ns1 = feature_ns1.sub_features[0]
        assert subfeature_ns1.type == "CDS"
        assert subfeature_ns1.qualifiers["gene"] == ["NS1"]
        assert subfeature_ns1.qualifiers["product"] == ["nonstructural protein 1"]
        assert subfeature_ns1.location.start == 24
        assert subfeature_ns1.location.end == 717

        assert feature_ns2.type == "gene"
        assert feature_ns2.qualifiers["gene"] == ["NS2"]
        # zero-based indexing for start using ExactPosition
        assert feature_ns2.location.start == 24
        assert feature_ns2.location.end == 862
        assert len(feature_ns2.sub_features) == 2
        subfeature_ns2_1 = feature_ns2.sub_features[0]
        subfeature_ns2_2 = feature_ns2.sub_features[1]
        assert subfeature_ns2_1.type == "CDS"
        assert subfeature_ns2_1.qualifiers["gene"] == ["NS2"]
        assert subfeature_ns2_1.qualifiers["product"] == ["nonstructural protein 2"]
        # check that subfeature location is correct
        assert subfeature_ns2_1.location.start == 24
        assert subfeature_ns2_1.location.end == 54

        assert subfeature_ns2_2.type == "CDS"
        assert subfeature_ns2_2.qualifiers["gene"] == ["NS2"]
        assert subfeature_ns2_2.qualifiers["product"] == ["nonstructural protein 2"]
        # check that subfeature location is correct
        assert subfeature_ns2_2.location.start == 526
        assert subfeature_ns2_2.location.end == 862
