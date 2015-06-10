

    import os
    
    !rm -r book
    os.makedirs('./book')
    
    cmd = "ipython nbconvert %s --to markdown --output %s"
    
    for root, dirs, files in os.walk('./'):
        for f in files:
            if '.ipynb_checkpoints' in root: continue
            fn, ext = os.path.splitext(f)
            if ext =='.ipynb':
                try:
                    os.makedirs(os.path.join('book', root))
                except OSError:
                    pass
                inf = os.path.join(root, f)
                outf = os.path.join('book', root, '%s.md' % fn)
                x = cmd % (inf, outf)
                print x
                    

    ipython nbconvert ./Index.ipynb --to markdown --output book/./Index.md
    ipython nbconvert ./Untitled.ipynb --to markdown --output book/./Untitled.md
    ipython nbconvert ./applications/biological-diversity.ipynb --to markdown --output book/./applications/biological-diversity.md
    ipython nbconvert ./fundamentals/database-searching.ipynb --to markdown --output book/./fundamentals/database-searching.md
    ipython nbconvert ./fundamentals/msa-assignment.ipynb --to markdown --output book/./fundamentals/msa-assignment.md
    ipython nbconvert ./fundamentals/multiple-sequence-alignment.ipynb --to markdown --output book/./fundamentals/multiple-sequence-alignment.md
    ipython nbconvert ./fundamentals/pairwise-alignment-exercises.ipynb --to markdown --output book/./fundamentals/pairwise-alignment-exercises.md
    ipython nbconvert ./fundamentals/pairwise-alignment.ipynb --to markdown --output book/./fundamentals/pairwise-alignment.md
    ipython nbconvert ./fundamentals/phylogeny-reconstruction.ipynb --to markdown --output book/./fundamentals/phylogeny-reconstruction.md
    ipython nbconvert ./fundamentals/sequence-mapping-and-clustering.ipynb --to markdown --output book/./fundamentals/sequence-mapping-and-clustering.md
    ipython nbconvert ./getting-started/reading-iab.ipynb --to markdown --output book/./getting-started/reading-iab.md



    !rm -r book


    
