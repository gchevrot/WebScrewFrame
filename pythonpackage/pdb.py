import nglview
from biopandas.pdb import PandasPdb
from .read import read_pdbid
from .data import aa_sidechain_chemical_properties as aa_hydrophobicity

class Pdb:
    """
    From a pdb id, the ScrewFrame parameters are computed
    """
    def __init__(self, pdbid):
        if pdbid[-4:] == ".pdb":
            self.pdbid = pdbid[:-4]
        else:
            self.pdbid = pdbid
        # raw pdb data
        self.pdb = read_pdbid(pdbid)
        # pandas pdb (biopandas)
        self.ppdb = PandasPdb().fetch_pdb(self.pdbid)
        self._coordinates = None

    def show(self):
        """return a 3D view of the protein within the notebook"""
        view = nglview.show_pdbid(self.pdbid)
        return view

    def atoms(self):
        """
        Extracting the lines starting with ATOM
        If there are 2 locations for the same atoms (column alt_loc is not empty),
        only the first location with alt_loc=='A' is used

        :return: lines starting with ATOM
        :rtype: pd.DataFrame
        """
        atoms = self.ppdb.df['ATOM']
        return atoms[(atoms['alt_loc'] == '') | (atoms['alt_loc'] == 'A')]

    def ca(self, chain_id_as_index=True):
        """
        Extracting carbon alpha ca_coordinates

        :param chain_id_as_index: if True the DataFrame indexes correspond to
                                  the chain_id
        :type chain_id_as_index: bool
        :return: carbon alpha coordinates
        :rtype: pd.DataFrame
        """
        df_atoms = self.atoms()
        if chain_id_as_index:
            self._coordinates = df_atoms[df_atoms['atom_name'] == 'CA'].filter(
                            ['chain_id', 'x_coord', 'y_coord', 'z_coord'],
                            axis='columns').set_index('chain_id')
        else:
            self._coordinates = df_atoms[df_atoms['atom_name'] == 'CA'].filter(
                            regex='coord$', axis='columns')
        return self

    def coordinates(self):
        """
        Return the coordinates defined by another method (e.g. the ca() method
        will return the carbon-alpha coordinates)
        Chained method to ca() method
        """
        return self._coordinates

    def distances(self):
        """
        Distances between the successive atoms
        Chained method to ca() method
        """
        # Distances between atoms within a same chain_id
        if self._coordinates.index.name == 'chain_id':
            return self._coordinates.groupby(by='chain_id').diff().dropna().pow(2).sum(axis='columns').pow(1/2)
        # Distances between successive atoms
        else:
            return self._coordinates.diff().dropna().pow(2).sum(axis='columns').pow(1/2)

    def description(self):
        """
        Number of chains and residues (in each chain)
        """
        df_atoms = self.atoms()
        # with groupby, counting the number of 'x_coord' for each residue
        df = df_atoms[df_atoms['atom_name'] == 'CA'].filter(regex='coord$',
                        axis='columns').groupby(by=df_atoms['chain_id']).count()
        n_chains = len(df.index)
        chain_names = [name for name in df.index]
        chain_lengths = [str(length) for length in df['x_coord']]
        print(f"""The molecule contains {n_chains} chains ({', '.join(chain_names)}) with {', '.join(chain_lengths)} residues respectively.""")

    def hydrophobicity(self):
        """
        Return the percentage of hydrophobic residue and the dataframe with the
        number of hydrophobic / hydrophilic residues

        The number of residue is evaluate with the CA atoms
        We map the hydrophobicity in function of the residue name
        and compute the percentage of hydrophobic residues
        """
        residues = self.atoms()[self.atoms()['atom_name'] == 'CA']['residue_name']
        # map the hydrophobicity and
        # count the number of hydrophbic and hydrophilic residues
        residues = residues.map(aa_hydrophobicity).value_counts(dropna=True)
        return residues['hydrophobic'] / residues.sum() * 100, residues
