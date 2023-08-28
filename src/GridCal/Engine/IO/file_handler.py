# GridCal
# Copyright (C) 2015 - 2023 Santiago Peñate Vera
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 3 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
import os
import json

from typing import Callable, Union, List
from GridCal.Engine.basic_structures import Logger
from GridCal.Engine.data_logger import DataLogger

from GridCal.Engine.IO.gridcal.json_parser import save_json_file_v3
from GridCal.Engine.IO.cim.cim16.cim_parser import CIMExport
from GridCal.Engine.IO.gridcal.excel_interface import save_excel, load_from_xls, interpret_excel_v3, interprete_excel_v2
from GridCal.Engine.IO.gridcal.pack_unpack import create_data_frames, data_frames_to_circuit
from GridCal.Engine.IO.matpower.matpower_parser import interpret_data_v1
from GridCal.Engine.IO.power_factory.dgs_parser import dgs_to_circuit
from GridCal.Engine.IO.matpower.matpower_parser import parse_matpower_file
from GridCal.Engine.IO.others.dpx_parser import load_dpx
from GridCal.Engine.IO.others.ipa_parser import load_iPA
from GridCal.Engine.IO.gridcal.json_parser import parse_json, parse_json_data_v2, parse_json_data_v3
from GridCal.Engine.IO.psse.raw_parser import read_raw
from GridCal.Engine.IO.psse.psse_to_gridcal import psse_to_gridcal
from GridCal.Engine.IO.power_world.power_world_parser import PowerWorldParser
# from GridCal.Engine.IO.cim.cim16.cim_parser import CIMImport
from GridCal.Engine.IO.cim.cgmes_2_4_15.cgmes_circuit import CgmesCircuit
from GridCal.Engine.IO.cim.cgmes_2_4_15.cgmes_to_gridcal import cgmes_to_gridcal
from GridCal.Engine.IO.gridcal.zip_interface import save_data_frames_to_zip, get_frames_from_zip
from GridCal.Engine.IO.gridcal.sqlite_interface import save_data_frames_to_sqlite, open_data_frames_from_sqlite
from GridCal.Engine.IO.gridcal.h5_interface import save_h5, open_h5
from GridCal.Engine.IO.psse.rawx_parser import rawx_parse, rawx_writer
from GridCal.Engine.IO.pypsa.pypsa_parser import parse_netcdf, parse_hdf5
from GridCal.Engine.Core.Devices.multi_circuit import MultiCircuit


class FileOpen:
    """
    File open interface
    """

    def __init__(self, file_name: Union[str, List[str]]):
        """
        File open handler
        :param file_name: name of the file
        """
        self.file_name: Union[str, List[str]] = file_name

        self.circuit: Union[MultiCircuit, None] = None

        self.cgmes_circuit: Union[CgmesCircuit, None] = None

        self.json_files = dict()

        self.logger = Logger()

        self.cgmes_logger = DataLogger()

    def open(self, text_func: Union[None, Callable] = None, progress_func: Union[None, Callable] = None):
        """
        Load GridCal compatible file
        :param text_func: pointer to function that prints the names
        :param progress_func: pointer to function that prints the progress 0~100
        :return: logger with information
        """
        self.logger = Logger()

        if isinstance(self.file_name, list):

            for f in self.file_name:
                _, file_extension = os.path.splitext(f)
                if file_extension.lower() not in ['.xml', '.zip']:
                    raise Exception('Loading multiple files that are not XML/Zip (xml or zip is for CIM or CGMES)')

            self.cgmes_circuit = CgmesCircuit(text_func=text_func, progress_func=progress_func, logger=self.cgmes_logger)
            self.cgmes_circuit.parse_files(cim_files=self.file_name)
            self.circuit = cgmes_to_gridcal(cgmes_model=self.cgmes_circuit, logger=self.cgmes_logger)

        else:

            if os.path.exists(self.file_name):
                name, file_extension = os.path.splitext(self.file_name)
                # print(name, file_extension)
                if file_extension.lower() in ['.xls', '.xlsx']:

                    data_dictionary = load_from_xls(self.file_name)

                    # Pass the table-like data dictionary to objects in this circuit
                    if 'version' not in data_dictionary:
                        interpret_data_v1(self.circuit, data_dictionary, self.logger)

                    elif data_dictionary['version'] == 2.0:
                        interprete_excel_v2(self.circuit, data_dictionary)

                    elif data_dictionary['version'] == 3.0:
                        interpret_excel_v3(self.circuit, data_dictionary)

                    elif data_dictionary['version'] == 4.0:
                        if data_dictionary is not None:
                            self.circuit = data_frames_to_circuit(data_dictionary, logger=self.logger)
                        else:
                            self.logger.add("Error while reading the file :(")
                            return None

                    else:
                        self.logger.add('The file could not be processed')

                elif file_extension.lower() == '.gridcal':

                    # open file content
                    data_dictionary, self.json_files = get_frames_from_zip(self.file_name,
                                                                           text_func=text_func,
                                                                           progress_func=progress_func,
                                                                           logger=self.logger)
                    # interpret file content
                    if data_dictionary is not None:
                        self.circuit = data_frames_to_circuit(data_dictionary, logger=self.logger)
                    else:
                        self.logger.add("Error while reading the file :(")
                        return None

                elif file_extension.lower() == '.sqlite':

                    # open file content
                    data_dictionary = open_data_frames_from_sqlite(self.file_name,
                                                                   text_func=text_func,
                                                                   progress_func=progress_func)
                    # interpret file content
                    if data_dictionary is not None:
                        self.circuit = data_frames_to_circuit(data_dictionary, logger=self.logger)
                    else:
                        self.logger.add("Error while reading the file :(")
                        return None

                elif file_extension.lower() == '.dgs':
                    self.circuit = dgs_to_circuit(self.file_name)

                elif file_extension.lower() == '.gch5':
                    self.circuit = open_h5(self.file_name, text_func=text_func, prog_func=progress_func)

                elif file_extension.lower() == '.m':
                    self.circuit, log = parse_matpower_file(self.file_name)
                    self.logger += log

                elif file_extension.lower() == '.dpx':
                    self.circuit, log = load_dpx(self.file_name)
                    self.logger += log

                elif file_extension.lower() == '.json':

                    # the json file can be the GridCal one or the iPA one...
                    with open(self.file_name, encoding="utf-8") as f:
                        data = json.load(f)

                        if isinstance(data, dict):
                            if 'Red' in data.keys():
                                self.circuit = load_iPA(self.file_name)
                            elif sum([x in data.keys() for x in ['type', 'version']]) == 2:
                                version = int(float(data['version']))

                                if data['type'] == 'Grid Exchange Json File' and 'profiles' in data.keys():
                                    if version == 2:
                                        self.circuit = parse_json_data_v2(data, self.logger)

                                    elif version == 3:
                                        self.circuit = parse_json_data_v3(data, self.logger)
                                    else:
                                        self.logger.add_error('Recognised as a gridCal compatible Json '
                                                              'but the version is not supported')

                            else:
                                self.logger.add_error('Unknown json format')

                        elif isinstance(data, list):
                            self.circuit = parse_json(self.file_name)

                        else:
                            self.logger.add_error('Unknown json format')

                elif file_extension.lower() == '.ejson3':
                    with open(self.file_name, encoding="utf-8") as f:
                        data = json.load(f)
                        self.circuit = parse_json_data_v3(data, self.logger)

                elif file_extension.lower() == '.raw':
                    pss_grid = read_raw(self.file_name,
                                        text_func=text_func,
                                        progress_func=progress_func,
                                        logger=self.logger)
                    self.circuit = psse_to_gridcal(psse_circuit=pss_grid, logger=self.logger)

                elif file_extension.lower() == '.rawx':
                    pss_grid = rawx_parse(self.file_name, logger=self.logger)
                    self.circuit = psse_to_gridcal(psse_circuit=pss_grid, logger=self.logger)

                elif file_extension.lower() == '.epc':
                    parser = PowerWorldParser(self.file_name)
                    self.circuit = parser.circuit
                    self.logger += parser.logger

                elif file_extension.lower() in ['.xml', '.zip']:
                    self.cgmes_circuit = CgmesCircuit(text_func=text_func,
                                                      progress_func=progress_func,
                                                      logger=self.cgmes_logger)
                    self.cgmes_circuit.parse_files(cim_files=[self.file_name])
                    self.circuit = cgmes_to_gridcal(cgmes_model=self.cgmes_circuit, logger=self.cgmes_logger)

                elif file_extension.lower() == '.hdf5':
                    self.circuit = parse_hdf5(self.file_name, self.logger)

                elif file_extension.lower() == '.nc':
                    self.circuit = parse_netcdf(self.file_name, self.logger)

            else:
                # warn('The file does not exist.')
                self.logger.add_error('Does not exist', self.file_name)

        return self.circuit

    def check_json_type(self, file_name):
        """
        Check the json file type from its internal data
        :param file_name: file path
        :return: data['type'] | 'Not json file' | 'Unknown json file'
        """

        if not os.path.exists(file_name):
            return 'Not json file'

        _, file_extension = os.path.splitext(file_name)

        if file_extension.lower() == '.json':
            with open(self.file_name, encoding="utf-8") as f:
                data = json.load(f)
                return data['type'] if 'type' in data else 'Unknown json file'
        else:
            return 'Not json file'


class FileSave:
    """
    FileSave
    """

    def __init__(self, circuit: MultiCircuit, file_name, text_func=None, progress_func=None,
                 simulation_drivers=list(), sessions=list(), json_files=dict()):
        """
        File saver
        :param circuit: MultiCircuit
        :param file_name: file name to save to
        :param text_func: Pointer to the text function
        :param progress_func: Pointer to the progress function
        :param simulation_drivers: List of Simulation Drivers
        :param sessions: List of sessions
        :param json_files: Dictionary of json files
        """
        self.circuit = circuit

        self.file_name = file_name

        self.simulation_drivers = simulation_drivers

        self.sessions = sessions

        self.json_files = json_files

        self.text_func = text_func

        self.progress_func = progress_func

    def save(self):
        """
        Save the file in the corresponding format
        :return: logger with information
        """
        if self.file_name.endswith('.xlsx'):
            logger = self.save_excel()

        elif self.file_name.endswith('.gridcal'):
            logger = self.save_zip()

        elif self.file_name.endswith('.sqlite'):
            logger = self.save_sqlite()

        elif self.file_name.endswith('.ejson3'):
            logger = self.save_json_v3()

        elif self.file_name.endswith('.xml'):
            logger = self.save_cim()

        elif self.file_name.endswith('.gch5'):
            logger = self.save_h5()

        elif self.file_name.endswith('.rawx'):
            logger = self.save_rawx()

        elif self.file_name.endswith('.newton'):
            logger = self.save_newton()

        elif self.file_name.endswith('.pgm'):
            logger = self.save_pgm()

        else:
            logger = Logger()
            logger.add_error('File path extension not understood', self.file_name)

        return logger

    def save_excel(self):
        """
        Save the circuit information in excel format
        :return: logger with information
        """

        logger = save_excel(self.circuit, self.file_name)

        return logger

    def save_zip(self):
        """
        Save the circuit information in zip format
        :return: logger with information
        """

        logger = Logger()

        dfs = create_data_frames(self.circuit)

        save_data_frames_to_zip(dfs=dfs,
                                filename_zip=self.file_name,
                                sessions=self.sessions,
                                diagrams=self.circuit.diagrams,
                                json_files=self.json_files,
                                text_func=self.text_func,
                                progress_func=self.progress_func)

        return logger

    def save_sqlite(self):
        """
        Save the circuit information in sqlite
        :return: logger with information
        """

        logger = Logger()

        dfs = create_data_frames(self.circuit)

        save_data_frames_to_sqlite(dfs,
                                   file_path=self.file_name,
                                   text_func=self.text_func,
                                   progress_func=self.progress_func)

        return logger

    def save_json_v3(self):
        """
        Save the circuit information in json format
        :return:logger with information
        """

        logger = save_json_file_v3(self.file_name, self.circuit, self.simulation_drivers)
        return logger

    def save_cim(self):
        """
        Save the circuit information in CIM format
        :return: logger with information
        """

        cim = CIMExport(self.circuit)
        cim.save(file_name=self.file_name)

        return cim.logger

    def save_h5(self):
        """
        Save the circuit information in CIM format
        :return: logger with information
        """

        logger = save_h5(self.circuit, self.file_name,
                         text_func=self.text_func,
                         prog_func=self.progress_func)

        return logger

    def save_rawx(self):
        """
        Save the circuit information in json format
        :return:logger with information
        """

        logger = rawx_writer(self.file_name, self.circuit)
        return logger

    def save_newton(self):
        """
        Save the circuit information in sqlite
        :return: logger with information
        """
        from GridCal.Engine.Core.Compilers.circuit_to_newton_pa import to_newton_pa, npa
        logger = Logger()

        time_series = self.circuit.time_profile is not None

        if time_series:
            tidx = list(range(len(self.circuit.time_profile)))
        else:
            tidx = None

        newton_grid, _ = to_newton_pa(self.circuit, time_series=time_series, time_indices=tidx)

        npa.FileHandler().save(newton_grid, self.file_name)

        return logger

    def save_pgm(self):
        """
        Save to Power Grid Model format
        :return: logger with information
        """
        from GridCal.Engine.Core.Compilers.circuit_to_pgm import save_pgm
        logger = Logger()

        save_pgm(filename=self.file_name, circuit=self.circuit, logger=logger, time_series=self.circuit.has_time_series)

        return logger
