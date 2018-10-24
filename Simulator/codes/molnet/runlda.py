LDA_CODE_PATH = '/Users/simon/git/lda/code'
MOTIF_DB_PATH = '/home/simon/git/motifdb'



import sys,csv,getopt,os


options = {
    'K': 20,
    'n_its':100,
    'overlap_thresh':0.3,
    'p_thresh':0.3,
    'lda_path':LDA_CODE_PATH,
    'pdf_backend':'Agg',
    'motifdb': '',
    'motifdb_path':MOTIF_DB_PATH,
}


def check_edge(edge,edges):
    if check_uni_edge(edge[0],edge[1],edges):
        return True
    if check_uni_edge(edge[1],edge[0],edges):
        return True
    return False

def check_uni_edge(node1,node2,edges):
    node1edges = [x for x in edges if x[0] == node1]
    if len(node1edges) > 0:
        node12edges = [x for x in node1edges if x[1] == node2]
        if len(node12edges) > 0:
            return True
    return False



if __name__ == '__main__':


    long_options = [
                'K=',
                'n_its=',
                'overlap_thresh=',
                'p_thresh=',
                'lda_path=',
                'pdf_backend=',
                'motifdb=',
                'motifdb_path=',
    ]

    

    input_prefix = sys.argv[1]

    found_options,the_rest = getopt.getopt(sys.argv[2:],"",long_options)


    for key,value in found_options:
        keyk = key.split('--')[1]
        try:
            options[keyk] = float(value)
        except:
            options[keyk] = value

    print(options)

    sys.path.append(options['lda_path'])

    from lda import VariationalLDA,compute_overlap_scores,write_summary_file,write_topic_report
    from ms2lda_feature_extraction import LoadMGF,MakeBinnedFeatures

    mgf_file = input_prefix + '.mgf'

    l = LoadMGF()
    
    print("Loading data from {}".format(mgf_file))
    
    ms1,ms2,metadata = l.load_spectra([mgf_file])

    print("Loaded {} spectra".format(len(ms1)))

    print("Creating corpus")

    corpus,word_mz_range = MakeBinnedFeatures().make_features(ms2)



    if len(options['motifdb'])>0: 
        print("Loading motifs from motifdb")

        db_path = options['motifdb_path'] + os.sep + 'motifs'
        db_code = options['motifdb_path'] + os.sep + 'code' + os.sep + 'utilities'   
        
        sys.path.append(db_code)

        from motifdb_loader import load_db,MotifFilter
        
        db_list = options['motifdb'].split()
        db_spectra,db_metadata = load_db(db_list,db_path)
        db_spectra,db_metadata = MotifFilter(db_spectra,db_metadata).filter()
    else:
        db_spectra = None
        db_metadata = None




    print("Creating LDA object with K = {}".format(options['K']))
    v = VariationalLDA(corpus[list(corpus.keys())[0]],
        K = int(options['K']),
        normalise = 1000.0,
        fixed_topics = db_spectra,
        fixed_topics_metadata = db_metadata)

    print("Running LDA for {} iterations".format(options['n_its']))

    v.run_vb(initialise = True,n_its = int(options['n_its']))

    print("Writing dictionary")
    output_dict_file = input_prefix + '_lda.dict'
    vd = v.make_dictionary(metadata = metadata,
        features = word_mz_range,
        compute_overlaps = True,
        filename = output_dict_file)

    print("Augmenting edge file")

    edge_file = input_prefix + '_edges.csv'
    original_edges = []
    with open(edge_file,'r') as f:
        reader = csv.reader(f)
        for edge in reader:
            original_edges.append(edge[:-1] + ['cosine'] + edge[-1:])

    print("Computing MS2LDA edges with p_thresh = {}, and overlap_thresh = {}".format(
        options['p_thresh'],options['overlap_thresh']))

    motifs = list(vd['beta'].keys())
    docs = list(vd['theta'].keys())
    motif_docs = {m: [] for m in motifs}
    for m in motifs:
        for doc in docs:
            if vd['theta'][doc].get(m,0.0) >= options['p_thresh'] and vd['overlap_scores'][doc].get(m,0.0) >= options['overlap_thresh']:
                motif_docs[m].append(doc)

    edges = []
    for m,docs in list(motif_docs.items()):
        for i,doc in enumerate(docs[:-1]):
            for doc2 in docs[i+1:]:
                edges.append([metadata[doc]['cid'],metadata[doc2]['cid'],m,0.0])


    # Write a file with all of the edges, even the ones across motifs
    all_edges = original_edges + edges
    with open(input_prefix+'_all_edges_ms2lda.csv','w') as f:
        writer = csv.writer(f)
        for line in all_edges:
            writer.writerow(line)


    # Filter out motif edges that do not have a cosine edge -- i.e. only keep those in the same mgf_file
    
    new_edges = []
    for edge in edges:
        if check_edge(edge,original_edges):
            new_edges.append(edge)
    edges = new_edges

    print("Combining edges")

    all_edges = original_edges + edges

    print("Writing combined edges")

    with open(input_prefix+'_edges_ms2lda.csv','w') as f:
        writer = csv.writer(f)
        for line in all_edges:
            writer.writerow(line)

    print("Creating summary file")

    summary_file = input_prefix + '_lda_summary.csv'
    write_summary_file(vd,summary_file)


    print("Creating pdf topic report")
    report_file = input_prefix + '_lda_report.pdf'
    write_topic_report(vd,report_file,backend = 'Agg')
            
def check_edge(edge,edges):
    if check_uni_edge(edge[0],edge[1],edges):
        return True
    if check_uni_edge(edge[1],edge[0],edges):
        return True
    return False

def check_uni_edge(node1,node2,edges):
    node1edges = list(filter(lambda x: x[0] == node1))
    if len(node1edges) > 0:
        node12edges = list(filter(lambda x: x[1] == node2))
        if len(node12edges) > 0:
            return True
    return False

