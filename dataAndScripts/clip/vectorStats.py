import pandas

df = pandas.read_excel('irCLIP_log.xlsx')
def n_genes(_df):
	genes = set([x.strip() for x in _df['Gene'].tolist()])
	print """\t{0} genes, {1} constructs.""".format(
		len(genes), len(_df.index))
	print ", ".join(genes)
print "All vectors:"
n_genes(df)
for vect in set(df['Vector'].tolist()):
	subdf = df[df['Vector']==vect].copy()
	print "Vector: {0}".format(vect)
	n_genes(subdf)
	print '-'

df = pandas.read_excel('irCLIP_log.xlsx', sheetname='Infections')
print 'Infections:'
n_genes(df)


