import gtf




def test_none_line():
    assert gtf.parse(None) is None


def test_comment_line():
    assert gtf.parse("# some comment") is None
    assert gtf.parse("  # some comment") is None


def test_feature_inclusion():
    line = """chr1\tHAVANA\tgene\t11869\t14412\t.\t+\t."""
    f = gtf.parse(line, features=["gene"])
    assert f is not None


def test_feature_exclusion():
    line = """chr1\tHAVANA\tgene\t11869\t14412\t.\t+\t."""
    f = gtf.parse(line, features=["exon"])
    assert f is None


def test_line_parsing():
    line = """chr1\tHAVANA\tgene\t11869\t14412\t.\t+\t."""
    f = gtf.parse(line)
    assert f is not None
    assert f.seqname == "chr1"
    assert f.source == "HAVANA"
    assert f.feature == "gene"
    assert f.start == 11869
    assert f.end == 14412
    assert f.score == -1
    assert f.strand == 1
    assert f.frame == -1


def test_line_parsing_attributes():
    line = """chr1\tHAVANA\tgene\t11869\t14412\t.\t+\t.\tgene_id "ENSG00000223972.4"; transcript_id "ENSG00000223972.4"; gene_type "pseudogene"; gene_status "KNOWN"; gene_name "DDX11L1"; transcript_type "pseudogene"; transcript_status "KNOWN"; transcript_name "DDX11L1"; level 2; havana_gene "OTTHUMG00000000961.2";"""
    f = gtf.parse(line)
    assert f.attributes["gene_id"] == "ENSG00000223972.4"
    assert f.attributes["transcript_id"] == "ENSG00000223972.4"
    assert f.attributes["gene_type"] == "pseudogene"
    assert f.attributes["gene_status"] == "KNOWN"
    assert f.attributes["gene_name"] == "DDX11L1"
    assert f.attributes["transcript_type"] == "pseudogene"
    assert f.attributes["transcript_status"] == "KNOWN"
    assert f.attributes["transcript_name"] == "DDX11L1"
    assert f.attributes["level"] == 2
    assert f.attributes["havana_gene"] == "OTTHUMG00000000961.2"


def test_line_parsing_attributes_include():
    line = """chr1\tHAVANA\tgene\t11869\t14412\t.\t+\t.\tgene_id "ENSG00000223972.4"; transcript_id "ENSG00000223972.4"; gene_type "pseudogene"; gene_status "KNOWN"; gene_name "DDX11L1"; transcript_type "pseudogene"; transcript_status "KNOWN"; transcript_name "DDX11L1"; level 2; havana_gene "OTTHUMG00000000961.2";"""
    f = gtf.parse(line, attributes=["gene_id"])
    assert f.attributes["gene_id"] == "ENSG00000223972.4"
    assert len(f.attributes) == 1

