import pandas as pd
import maboss

class Cellline:
    def __init__(self, name, mutations, transition_rates_up, initial_states, dict_gene_nodes, dict_strict_gene_nodes):
        self.name = name
        self.mutations = mutations
        self.transition_rates_up = transition_rates_up
        self.initial_states = initial_states
        self.dict_gene_nodes = dict_gene_nodes
        self.dict_strict_gene_nodes = dict_strict_gene_nodes
    def __repr__(self):
        try :
            self.model
            try:
                self.results
                return self.name+', simu'
            except:
                return self.name+', compu'
        except:
            return self.name 
    
    def personalize_model(self, model):
        personalized_model = model.copy()
        nodes = model.network.names
        
        for node in nodes:
            to_inverse=False
            if node in self.dict_gene_nodes:
                for gene in self.dict_gene_nodes[node]:
                    if gene[0]=='!':
                        to_inverse = True
                        gene = gene[1:]
                    if gene in self.mutations.keys():
                        val_mut = self.mutations[gene]
                        if to_inverse:
                            if val_mut=='ON':
                                val_mut ='OFF'
                            else:
                                val_mut = 'ON'
                        personalized_model.mutate(node, val_mut)
            if node in self.dict_strict_gene_nodes:
                for gene in self.dict_strict_gene_nodes[node]:
                    if gene in self.transition_rates_up.keys():
                        personalized_model.param['$u_'+node] = self.transition_rates_up[gene]
                        personalized_model.param['$d_'+node] = 1/self.transition_rates_up[gene]
            if len(self.initial_states)!=0:
                for node in self.initial_states:
                    personalized_model.network.set_istate(node, self.initial_states[node])
        self.model = personalized_model
        return personalized_model
    
    def run_simulation(self):
        self.results = self.model.run()
        print('done')

class CellEnsemble:
    def __init__(self, celllines, general_model):
        self.model = general_model
        self.data = celllines
        self.simu = {'base':{cellname:Cellline(**celllines[cellname])
                             for cellname in celllines}}
        self.update()
    def __repr__(self):
        return 'conditions : '+str(list(self.simu.keys()))+'\ncell lines : '+str(list(self.simu['base'].keys()))
    def resume(self):
        """Plot the names of the personalized models stored, under a pandas.DataFrame format."""
        return pd.DataFrame(self.simu)
    
    def add_condition(self, conditions, condition_name):
        """Add a condition (new version of each cell line personalized model with mutations) to the CellEnsembl object
        conditions : mutations in the MaBoSS format : (Node_name, effect). The effect can be 'ON' or 'OFF'. For mutliple mutation, a list or tuple can be given.
        condition_name : name of the condition to call them and to define title of the output graphs.
        """
        self.simu[condition_name] = {cellname:Cellline(**celllines[cellname])
                                     for cellname in celllines}
        self.update()
        if type(conditions[0]) is not str:
            print('{} mutations for the condition {}'.format(len(conditions), condition_name))
            for condition in conditions:
                for k in self.simu[condition_name]:
                    self.simu[condition_name][k].model.mutate(*condition)
        else:
            for k in self.simu[condition_name]:
                self.simu[condition_name][k].model.mutate(*conditions)
    
    def update(self, new_model=None):
        """By default, the method recreates personalized versions of the model which is stored as an attribute : self.model.
        new model : You can also give the new version of the general model directly as an argument in this method with this argument.
        It allow to keep the format of the CellEnsembl object (conditions and celllines) but to change the general model used.
        """
        if new_model is not None:
            self.model = new_model
            
        for k1 in self.simu:
            for k2 in self.simu[k1]:
                self.simu[k1][k2].personalize_model(self.model)
    
    def run_simulation(self, celllines='all', conditions='all', redo=True, mute=False):
        
        if celllines == 'all':
            celllines = list(self.simu['base'].keys())
        elif type(celllines) not in [list, tuple]:
            celllines = [celllines]
        if conditions == 'all':
            conditions = list(self.simu.keys())
        elif type(conditions) not in [list, tuple]:
            conditions = [conditions]
            
        for condition in conditions:
            for cellline in celllines:
                if redo==False:
                    try:
                        self.simu[condition][cellline].results
                        if mute==False:
                            print('the simulation {}|{} will not be re-computed'.format(cellline, condition))
                        continue
                    except:
                        pass
                if mute==False:
                    print('Simulating', cellline, 'in the condition', condition)
                self.simu[condition][cellline].run_simulation()
                
       
    def plot_piechart(self, celllines='all', conditions='all'):
        if celllines == 'all':
            celllines = list(self.simu['base'].keys())
        elif type(celllines) not in [list, tuple]:
            celllines = [celllines]
        if conditions == 'all':
            conditions = list(self.simu.keys())
        elif type(conditions) not in [list, tuple]:
            conditions = [conditions]
            
        
        if len(conditions)>1:
            for cellline in celllines:
                fig, axs = plt.subplots(1, len(conditions), figsize=(8*len(conditions),5), constrained_layout=True)
                for condition in conditions:
                    try: 
                        self.simu[condition][cellline].results.plot_piechart(axes = axs[conditions.index(condition)])
                        axs[conditions.index(condition)].set_title(cellline + ' | ' + condition)
                    except:
                        print('the simulation {}|{} seems not having been computed'.format(cellline, condition))
        else:
            for condition in conditions:
                for cellline in celllines:
                    try: 
                        self.simu[condition][cellline].results.plot_piechart()
                        plt.title(cellline + ' | ' + condition)
                    except:
                        print('the simulation {}|{} seems not having been computed'.format(cellline, condition))
    def plot_trajectory(self, celllines='all', conditions='all'):
        if celllines == 'all':
            celllines = list(self.simu['base'].keys())
        elif type(celllines) not in [list, tuple]:
            celllines = [celllines]
        if conditions == 'all':
            conditions = list(self.simu.keys())
        elif type(conditions) not in [list, tuple]:
            conditions = [conditions]
            
        
        if len(conditions)>1:
            for cellline in celllines:
                fig, axs = plt.subplots(1, len(conditions), figsize=(8*len(conditions),5), constrained_layout=True)
                for condition in conditions:
                    try: 
                        self.simu[condition][cellline].results.plot_trajectory(axes = axs[conditions.index(condition)])
                        axs[conditions.index(condition)].set_title(cellline + ' | ' + condition)
                    except:
                        print('the simulation {}|{} seems not having been computed'.format(cellline, condition))
        else:
            for condition in conditions:
                for cellline in celllines:
                    try: 
                        self.simu[condition][cellline].results.plot_trajectory()
                        plt.title(cellline + ' | ' + condition)
                    except:
                        print('the simulation {}|{} seems not having been computed'.format(cellline, condition))

                        
                        
    def compare_survival(self, celllines='all', conditions='all'):
        
        outputs_df = pd.DataFrame(index=[], columns=['proliferation_rates', 'survival_rates'])
        
        if celllines == 'all':
            celllines = list(self.simu['base'].keys())
        elif type(celllines) not in [list, tuple]:
            celllines = [celllines]
        if conditions == 'all':
            conditions = list(self.simu.keys())
        elif type(conditions) not in [list, tuple]:
            conditions = [conditions]
            
        for condition in conditions:
            for cellline in celllines:
                outputs_dict = {'survival_rates' : {}, 'proliferation_rates' : {}}

                st = self.simu[condition][cellline].results.last_states_probtraj
                death_rate = float(sum(st[s] for s in st if 'Apoptosis' in s or 'Mitotic' in s))
                proliferation = float(sum(st[s] for s in st if 'Apoptosis' not in s and 'Mitotic' not in s and 'Prolif' in s))

                outputs_dict['survival_rates'][condition+'_'+cellline] =  1-death_rate
                outputs_dict['proliferation_rates'][condition+'_'+cellline] = proliferation
                outputs_df = pd.concat([pd.DataFrame(outputs_dict), outputs_df])
        return outputs_df
                
                
