import py3Dmol
from math import ceil
import configparser
from IPython.display import HTML
from Bio.PDB import PDBParser, DSSP
import numpy as np

class CleavageSiteVisualizer:
    def __init__(self, config_path='config.ini'):
        """
        Args:
            config_path: Configuration file
        """
        self.config = configparser.ConfigParser()
        self.config.read(config_path)
        
    def _get_residue_coords(self, pdb_data, residue_num):
        for line in pdb_data.split('\n'):
            if line.startswith('ATOM') and int(line[22:26].strip()) == residue_num:
                return {
                    'x': float(line[30:38]),
                    'y': float(line[38:46]),
                    'z': float(line[46:54])
                }
        return None

    def _add_region_label(self, view, coords, label_text, background_color):
        view.addLabel(label_text, {
            'position': coords,
            'backgroundColor': background_color,
            'fontColor': self.config.get('Label_Colors', 'font_color'),
            'inFront': True,
            'fontSize': self.config.getint('Labels', 'font_size')
        })

    def _add_region_style(self, view, resi_range, style_type, color, opacity=None, scale=None):
        style_dict = {'color': color}
        if opacity is not None:
            style_dict['opacity'] = opacity
        if scale is not None:
            style_dict['scale'] = scale
            
        view.addStyle({'resi': resi_range}, {style_type: style_dict})

    def _calculate_rsa(self, pdb_file, rsa_cal):
        """Calculate RSA for all residues in the protein using Wilke values."""
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('protein', pdb_file)
        dssp = DSSP(structure[0], pdb_file, acc_array=rsa_cal)
        rsa_dict = {}
        for key in dssp.keys():
            residue_num = key[1][1]
            rsa = dssp[key][3]
            rsa_dict[residue_num] = rsa
        #print(rsa_dict)        
        return rsa_dict

    def _get_rsa_color(self, rsa):
        """Get color based on RSA value using a continuous blue-to-red gradient."""
        rsa = max(0, min(1, rsa))
        
        blue = np.array([0, 0, 255])
        red = np.array([255, 0, 0])
        
        color = blue + (red - blue) * rsa
        
        # Convert to hex format
        return f'#{int(color[0]):02x}{int(color[1]):02x}{int(color[2]):02x}'

    def visualize(self, pdb_path, extracellular_regions=None, transmembrane_regions=None, 
                 cleavage_sites=None, cleavage_sites_sheddase=None, 
                 cleavage_sites_ADAM=None, cleavage_sites_BACE=None,
                 show_rsa=True,rsa_cal='Wilke'):
        """
        Args:
        - pdb_path: pdb file
        - extracellular_regions: List of tuples [(start1, end1), (start2, end2), ...]
        - transmembrane_regions: List of tuples [(start1, end1), (start2, end2), ...]
        - cleavage_sites: List of all possible cleavage sites
        - cleavage_sites_sheddase: List of sheddase cleavage sites
        - cleavage_sites_ADAM: List of ADAM cleavage sites
        - cleavage_sites_BACE: List of BACE cleavage sites
        - show_rsa: Boolean to toggle RSA visualization (default: True)
        - rsa_cal: Accessible surface area (ASA) from either Miller et al. (1987), 
                   Sander & Rost (1994), 
                   or Wilke: Tien et al. 2013, as string Sander/Wilke/Miller. Defaults to Wilke.
        """
        with open(pdb_path, 'r') as f:
            pdb_data = f.read()
        
        view = py3Dmol.view(width=self.config.getint('Display', 'width'),
                           height=self.config.getint('Display', 'height'))
        view.addModel(pdb_data, 'pdb')
        
        # Base style
        view.setStyle({'cartoon': {'color': self.config.get('Colors', 'default_protein')}})
        
        # Calculate RSA if needed
        rsa_dict = self._calculate_rsa(pdb_path, rsa_cal) if show_rsa else {}
        
        # Base style with RSA coloring if enabled
        if show_rsa:
            for residue_num, rsa_value in rsa_dict.items():
                color = self._get_rsa_color(rsa_value)
                view.addStyle({'resi': str(residue_num)}, 
                            {'cartoon': {'color': color, 'opacity': 0.7}})
        
        # Extracellular regions
        if extracellular_regions:
            # Add label only for the first region
            first_region = extracellular_regions[0]
            middle_residue = ceil((first_region[0] + first_region[1]) / 2)
            coords = self._get_residue_coords(pdb_data, middle_residue)
            if coords:
                coords['x'] += self.config.getfloat('Labels', 'label_offset_x')
                coords['y'] += self.config.getfloat('Labels', 'label_offset_y')
                self._add_region_label(view, coords,
                                     self.config.get('Labels', 'extracellular_text'),
                                     self.config.get('Label_Colors', 'background_extracellular'))
            
            for start, end in extracellular_regions:
                if not show_rsa:  # Only apply region coloring if RSA is disabled
                    self._add_region_style(view,
                                         f'{start}-{end}',
                                         self.config.get('Styling', 'extracellular_style'),
                                         self.config.get('Colors', 'extracellular'),
                                         self.config.getfloat('Styling', 'extracellular_opacity'))

        # Transmembrane regions
        if transmembrane_regions:
            # Add label only for the first region
            first_region = transmembrane_regions[0]
            middle_residue = ceil((first_region[0] + first_region[1]) / 2)
            coords = self._get_residue_coords(pdb_data, middle_residue)
            if coords:
                coords['x'] += self.config.getfloat('Labels', 'label_offset_x')
                coords['y'] += self.config.getfloat('Labels', 'label_offset_y')
                self._add_region_label(view, coords,
                                     self.config.get('Labels', 'transmembrane_text'),
                                     self.config.get('Label_Colors', 'background_transmembrane'))

            for start, end in transmembrane_regions:
                if not show_rsa:  # Only apply region coloring if RSA is disabled
                    self._add_region_style(view,
                                         f'{start}-{end}',
                                         self.config.get('Styling', 'transmembrane_style'),
                                         self.config.get('Colors', 'transmembrane'),
                                         self.config.getfloat('Styling', 'transmembrane_opacity'))

        # Cleavage sites
        if cleavage_sites:
            cleavage_sites = set(cleavage_sites)
            sheddase_sites = set(cleavage_sites_sheddase if cleavage_sites_sheddase else [])
            adam_sites = set(cleavage_sites_ADAM if cleavage_sites_ADAM else [])
            bace_sites = set(cleavage_sites_BACE if cleavage_sites_BACE else [])
    
            # Validate relationships
            if adam_sites.intersection(bace_sites):
                raise ValueError("ADAM and BACE sites should not overlap")
            if not (adam_sites.issubset(sheddase_sites) and bace_sites.issubset(sheddase_sites)):
                raise ValueError("ADAM and BACE sites must be subsets of sheddase sites")
            if not sheddase_sites.issubset(set(cleavage_sites)):
                raise ValueError("Sheddase sites must be a subset of cleavage sites")
    
            labeled_sites = set()
    
            # Add labels for each type (only once and at different sites)
            if adam_sites:
                # Find first ADAM site that isn't labeled yet
                for site in sorted(adam_sites):
                    if site not in labeled_sites:
                        coords = self._get_residue_coords(pdb_data, site)
                        if coords:
                            coords['x'] += self.config.getfloat('Labels', 'label_offset_x')
                            coords['y'] += self.config.getfloat('Labels', 'label_offset_y')
                            self._add_region_label(view, coords,
                                                 self.config.get('Labels', 'cleavage_text_ADAM'),
                                                 self.config.get('Label_Colors', 'background_cleavage_ADAM'))
                            labeled_sites.add(site)
                            break

            if bace_sites:
                # Find first BACE site that isn't labeled yet
                for site in sorted(bace_sites):
                    if site not in labeled_sites:
                        coords = self._get_residue_coords(pdb_data, site)
                        if coords:
                            coords['x'] += self.config.getfloat('Labels', 'label_offset_x')
                            coords['y'] += self.config.getfloat('Labels', 'label_offset_y')
                            self._add_region_label(view, coords,
                                                 self.config.get('Labels', 'cleavage_text_BACE'),
                                                 self.config.get('Label_Colors', 'background_cleavage_BACE'))
                            labeled_sites.add(site)
                            break

            if sheddase_sites - adam_sites - bace_sites:  # Only label sheddase sites that aren't ADAM or BACE
                # Find first sheddase site that isn't labeled yet
                for site in sorted(sheddase_sites - adam_sites - bace_sites):
                    if site not in labeled_sites:
                        coords = self._get_residue_coords(pdb_data, site)
                        if coords:
                            coords['x'] += self.config.getfloat('Labels', 'label_offset_x')
                            coords['y'] += self.config.getfloat('Labels', 'label_offset_y')
                            self._add_region_label(view, coords,
                                                 self.config.get('Labels', 'cleavage_text_sheddase'),
                                                 self.config.get('Label_Colors', 'background_cleavage_sheddase'))
                            labeled_sites.add(site)
                            break
                            
            if cleavage_sites - sheddase_sites:  # Only label cleavage sites that aren't sheddase_sites
                # Find first cleavage sites that isn't labeled yet
                for site in sorted(cleavage_sites - sheddase_sites):
                    if site not in labeled_sites:
                        coords = self._get_residue_coords(pdb_data, site)
                        if coords:
                            coords['x'] += self.config.getfloat('Labels', 'label_offset_x')
                            coords['y'] += self.config.getfloat('Labels', 'label_offset_y')
                            self._add_region_label(view, coords,
                                                 self.config.get('Labels', 'cleavage_text_other'),
                                                 self.config.get('Label_Colors', 'background_cleavage_other'))
                            labeled_sites.add(site)
                            break

            for site in cleavage_sites:
                if site in adam_sites:
                    style_color = self.config.get('Colors', 'cleavage_ADAM')
                elif site in bace_sites:
                    style_color = self.config.get('Colors', 'cleavage_BACE')
                elif site in sheddase_sites:
                    style_color = self.config.get('Colors', 'cleavage_sheddase')
                else:
                    style_color = self.config.get('Colors', 'cleavage_other')

                self._add_region_style(view,
                                     str(site),
                                     self.config.get('Styling', 'cleavage_style'),
                                     style_color,
                                     scale=self.config.getfloat('Styling', 'cleavage_sphere_scale'))

        view.zoomTo()
        view.render()
        return view

    def generate_html(self, pdb_path, extracellular_regions=None, transmembrane_regions=None, 
                     cleavage_sites=None, cleavage_sites_sheddase=None, 
                     cleavage_sites_ADAM=None, cleavage_sites_BACE=None,
                     show_rsa=True, div_id='pdb-viewer'):
        """
        Generate HTML code for protein visualization - **for future use in webapp**
        
        Args:
        [Previous args documentation...]
        - show_rsa: Boolean to toggle RSA visualization
        """
        view = self.visualize(pdb_path, 
                            extracellular_regions=extracellular_regions,
                            transmembrane_regions=transmembrane_regions,
                            cleavage_sites=cleavage_sites,
                            cleavage_sites_sheddase=cleavage_sites_sheddase,
                            cleavage_sites_ADAM=cleavage_sites_ADAM,
                            cleavage_sites_BACE=cleavage_sites_BACE,
                            show_rsa=show_rsa)

        return f"""
        <div id="{div_id}" style="height: {self.config.get('Display', 'height')}px; 
                                  width: {self.config.get('Display', 'width')}px; 
                                  position: relative;">
            <script>
                {view.js()}
            </script>
        </div>
        """