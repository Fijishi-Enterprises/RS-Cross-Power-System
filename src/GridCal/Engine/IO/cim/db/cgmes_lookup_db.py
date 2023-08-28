import os
from typing import List
from GridCal.Engine.IO.cim.db.base_db import BaseDb
from GridCal.Engine.IO.cim.db.file_system import get_create_roseta_db_folder
from GridCal.Engine.IO.cim.cgmes_2_4_15.cgmes_circuit import CgmesCircuit
from GridCal.Engine.IO.cim.cgmes_2_4_15.cim_parser import read_cgmes_files


class CgmesLookUpDb(BaseDb):
    """
    CgmesLookUpDb
    """

    def __init__(self, new_db=False):
        """

        :param new_db: By default, load from the disk
        """
        here = os.path.dirname(os.path.abspath(__file__))
        BaseDb.__init__(self,
                        db_folder=os.path.join(get_create_roseta_db_folder(), 'CgmesLookUp'),
                        db_extension='.zip',
                        new_db=new_db,
                        init_files=[os.path.join(here, '..', 'data', 'ENTSOe_boundary_set.zip')])

        self.circuit: CgmesCircuit | None = None
        print('loading CGMES boundary set...', end='')
        self.read_db_file(file_name=self.get_last_file_path())
        print('ok')

    def read_db_file(self, file_name):
        """
        Read CGMES database file
        :param file_name:
        :return:
        """
        pth = os.path.join(self.db_folder, file_name)
        self.circuit: CgmesCircuit = read_cgmes_files(cim_files=[pth])

    def get_structures_names(self) -> List[str]:
        classes = [prop.property_name for prop in self.circuit.get_class_properties()]
        classes.sort()

        return classes
