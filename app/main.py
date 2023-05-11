import pandas as pd


def scoring_algorithm(input_file, output_file, upper_threshold=.95, lower_threshold=.5):
    df_ls = []
    df = pd.read_csv(input_file, index_col=0).T
    for column in df:
        temp_df = df[column].dropna()
        loss_score = 0
        if len(temp_df) > 1:
            for index, row in temp_df.items():
                if row > upper_threshold:
                    loss_score += row
                elif row < lower_threshold:
                    loss_score -= 2
                else:
                    loss_score -= 1
            ls = [column, loss_score]
            df_ls.append(ls)
    pd.DataFrame(df_ls, columns=['gene_mutation', 'score']).to_csv(output_file, index=False)


def gene_enrichment(input_file, output_file):
    appearance_count, gene_score = {}, {}
    df = pd.read_csv(input_file, index_col=0)
    for i in df.index:
        gene = i.split('_')[0]
        if gene not in gene_score:
            gene_score[gene] = df.loc[i, ['Precision score', 'Specificity score']].sum()
            appearance_count[gene] = 1
        else:
            gene_score[gene] += df.loc[i, ['Precision score', 'Specificity score']].sum()
            appearance_count[gene] += 1
    updated_appearance_count = [key for key, value in appearance_count.items() if value > 1]
    df = pd.DataFrame().from_dict(gene_score, orient='index', columns=['Enriched gene score'])
    df.loc[[i for i in updated_appearance_count if i in df.index]].to_csv(output_file)
