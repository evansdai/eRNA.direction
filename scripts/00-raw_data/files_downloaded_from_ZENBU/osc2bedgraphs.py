import pandas as pd
from typing import Union

class osc2bedgraphs(object):
    def __init__(self, workfolder):
        self.expression_name = []
        self.workfolder = workfolder
        

    def read_osc_as_pd(self, filename:str, column_names_file: Union[str,None]):
        self.df = pd.read_csv(filename, sep='\t', )
        if not column_names_file is None:
            self.column_names = pd.read_csv(column_names_file, sep = '\t').columns
            self.df.columns = self.column_names
        self.df_plus = self.df[self.df['eedb:strand'] == '+']
        self.df_minus = self.df[self.df['eedb:strand'] == '-']
        assert(self.df[(self.df['eedb:strand'] != '+') & (self.df['eedb:strand'] != '-')].shape[0] ==0)
    
    def get_exp_name(self):
        colnames = self.df.columns
        for i in colnames:
            if i.startswith('exp.tagcount'):
                self.expression_name.append(i)
        # print(self.expression_name)

    def batch_save_bedgraph(self):
        for strand, df in [['+', self.df_plus],['-',self.df_minus]]:
            for i in self.expression_name:
                df_temp = df.loc[:, ['eedb:chrom','eedb:start.0base','eedb:end',i]]
                df_temp2 = df_temp.loc[df_temp[i]>0,:]
                print('{}.{}'.format(i, strand))
                df_temp2.to_csv('{}\{}.{}.bedGraph'.format(self.workfolder,i,strand), header = False, index = False, sep = '\t')


    # pool single cells
    def get_9_groups(self):
        df_summary = pd.read_csv(
            r"C:\Users\evans\Documents\albin_lab\single cell cage\data\sample_summary.csv", index_col='id')
        return {key: list(item.index) for key, item in df_summary.groupby('run_timepoint')}

    def get_3_groups(self):
        df_summary = pd.read_csv(
            r"C:\Users\evans\Documents\albin_lab\single cell cage\data\sample_summary.csv", index_col='id')
        return {key: list(item.index) for key, item in df_summary.groupby('timepoint')}

    def pool_cells_and_save_bedgraph(self):
        dict_groups = self.get_3_groups()
        for group in dict_groups.keys():
            ids = dict_groups[group]
            for strand, df in [['+', self.df_plus], ['-', self.df_minus]]:
                df_pooled_counts = df.loc[:, ['exp.tagcount_pm.{}'.format(id_i) for id_i in ids]].sum(axis = 1)
                df_pooled_counts.name = 'pooled_counts'
                df_seq = df.loc[:, ['eedb:chrom','eedb:start.0base','eedb:end']]
                df_temp = pd.concat([df_seq, df_pooled_counts], axis =1 )
                df_temp2 = df_temp.loc[df_temp['pooled_counts'] > 0, :]
                print('{}.{}'.format(group, strand))
                df_temp2.to_csv('{}\pooled_{}.{}.bedGraph'.format(self.workfolder,group,strand), header = False, index = False, sep = '\t')
    

if __name__ == '__main__':
    o2b = osc2bedgraphs(
        r"C:\Users\evans\Documents\albin_lab\single cell cage\Download raw reads\output\\")
    
    
    # Time course libraries normalised QCfiltered counts
    # filename = r"C:\Users\evans\Documents\albin_lab\single cell cage\data\bed\5BC1CAGE5DTGF-CEB2TimecourselibrariesnormalisedQCfilteredandgroupedbytimepoint_count.osc"
    # column_names_file = r"C:\Users\evans\Documents\albin_lab\single cell cage\data\bed\column_name2.osc"
    # o2b.read_osc_as_pd(filename, column_names_file)
    # o2b.get_exp_name()
    # o2b.batch_save_bedgraph()

    # Not used: Pool single cells by group
    # normalized and QC filtered counts
    # filename = r"C:\Users\evans\Documents\albin_lab\single cell cage\data\bed\5BC1CAGE5DTGF-CEB2TimecourselibrariesnormalisedQCfilteredandgroupedbytimepoint_count.osc"
    # column_names_file = r"C:\Users\evans\Documents\albin_lab\single cell cage\data\bed\column_name2.osc"
    # o2b.read_osc_as_pd(filename, column_names_file)
    # o2b.pool_cells_and_save_bedgraph()


    # Bulk data
    filename = r"C:\Users\evans\Documents\albin_lab\single cell cage\data\bed\5BBulkCAGE5DTGF-CEB2TimecoursebulklibrariesTPM-normalisedgroupedbytimepoint..osc"
    o2b.read_osc_as_pd(filename, column_names_file = None)
    o2b.get_exp_name()
    o2b.batch_save_bedgraph()
