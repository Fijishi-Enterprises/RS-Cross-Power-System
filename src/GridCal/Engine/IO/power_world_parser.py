# This file is part of GridCal.
#
# GridCal is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# GridCal is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GridCal.  If not, see <http://www.gnu.org/licenses/>.

import chardet
import os
from typing import List, AnyStr, Dict

import numpy as np

from GridCal.Engine.Core.multi_circuit import MultiCircuit
from GridCal.Engine.basic_structures import Logger
import GridCal.Engine.Devices as dev


def interpret_line(line, splitter=','):
    """
    Split text into arguments and parse each of them to an appropriate format (int, float or string)
    Args:
        line: text line
        splitter: value to split by
    Returns: list of arguments
    """
    parsed = list()
    elms = line.split(splitter)

    for elm in elms:
        try:
            # try int
            el = int(elm)
        except ValueError as ex1:
            try:
                # try float
                el = float(elm)
            except ValueError as ex2:
                # otherwise just leave it as string
                el = elm.strip()
        parsed.append(el)
    return parsed


def find_between(s, start, end):
    # if start in s and end in s:
    #     return (s.split(start))[1].split(end)[0]
    # else:
    #     return ""
    return (s.split(start))[1].split(end)[0]


def split_line(lne):

    chunks = list()

    if len(lne):

        current = ''
        started = False
        reading_str = False
        float_detected = False
        for chr in lne:
            if chr != ':':

                if chr != ' ' and not started:
                    # start any
                    reading_str = chr == '"'
                    started = True
                    float_detected = False
                    current += chr

                elif chr not in [' ', '"'] and started and not reading_str:
                    # keep reading value

                    if not float_detected:
                        if chr == '.':
                            float_detected = True

                    current += chr

                elif chr == ' ' and started and not reading_str:
                    # finalize reading value
                    started = False
                    if float_detected:
                        chunks.append(float(current))
                    else:
                        chunks.append(int(current))
                    current = ''

                elif chr != '"' and started and reading_str:
                    # keep reading string
                    current += chr

                elif chr == '"' and started and reading_str:
                    # finalize reading string
                    current += chr
                    started = False
                    chunks.append(current.replace('"', ''))
                    current = ''

        # add the last chunk
        if len(current):
            if reading_str:
                chunks.append(current)
            else:
                if float_detected:
                    chunks.append(float(current))
                else:
                    chunks.append(int(current))

    return chunks


def parse_substations(data_lst: List[List]):

    data = dict()
    for raw in data_lst:
        code = raw[0]
        data[code] = dev.Substation(name=raw[1],
                                    idtag=None,
                                    code=code,
                                    latitude=raw[2],
                                    longitude=raw[3])

    return data


def parse_buses(data_lst: List[List], substations_dict: Dict[int, dev.Substation]):
    buses_dict = dict()
    bus_volt = dict()
    for raw in data_lst:
        code = raw[0]
        area_idx = raw[9]
        zone_idx = raw[10]

        lat = raw[19]
        lon = raw[20]
        active = 1 - raw[18]

        if substations_dict is not None:
            st_idx = raw[26]  # or maybe 28
            substation = substations_dict[st_idx]

            if lat == 0:
                lat = substation.latitude

            if lon == 0:
                lon = substation.longitude
        else:
            substation = None

        bus_volt[code] = raw[6]

        buses_dict[code] = dev.Bus(name=raw[1],
                                   idtag=None,
                                   code=code,
                                   vnom=raw[2],
                                   vmin=0.9,
                                   vmax=1.1,
                                   angle_min=-6.28, angle_max=6.28, r_fault=0.0, x_fault=0.0,
                                   xpos=0, ypos=0, height=0, width=0,
                                   active=active,
                                   is_slack=False,
                                   is_dc=False,
                                   area=None,
                                   zone=None,
                                   substation=substation,
                                   country=None,
                                   longitude=lon,
                                   latitude=lat)

    return buses_dict, bus_volt


def parse_dc_buses(data_lst: List[List]):
    buses_dict = dict()
    bus_volt = dict()
    for raw in data_lst:
        code = raw[0]
        area_idx = raw[4]
        zone_idx = raw[5]

        lat = 0
        lon = 0
        active = True

        substation = None

        bus_volt[code] = raw[6] / raw[7]

        buses_dict[code] = dev.Bus(name=raw[1],
                                   idtag=None,
                                   code=code,
                                   vnom=raw[7],
                                   vmin=0.9,
                                   vmax=1.1,
                                   angle_min=-6.28, angle_max=6.28, r_fault=0.0, x_fault=0.0,
                                   xpos=0, ypos=0, height=0, width=0,
                                   active=active,
                                   is_slack=False,
                                   is_dc=True,
                                   area=None,
                                   zone=None,
                                   substation=substation,
                                   country=None,
                                   longitude=lon,
                                   latitude=lat)

    return buses_dict, bus_volt


def parse_transformers(data_lst: List[List], buses_dict: Dict[int, dev.Bus]):
    data = list()
    for raw in data_lst:
        name = '{0}_{1}_{2}'.format(raw[1], raw[4], raw[6])
        code = '{0}_{1}_{2}'.format(raw[0], raw[3], raw[6])
        bus_f = buses_dict[raw[0]]
        bus_t = buses_dict[raw[3]]
        Vh = max(raw[29], raw[30])
        Vl = min(raw[29], raw[30])
        rate = raw[35]
        r = raw[23]
        x = raw[24]
        elm = dev.Transformer2W(bus_from=bus_f,
                                bus_to=bus_t,
                                HV=Vh,
                                LV=Vl,
                                name=name,
                                idtag=None,
                                code=code,
                                r=r,
                                x=x,
                                g=1e-20,
                                b=1e-20,
                                rate=rate,
                                tap=1.0,
                                tap_module_max=1.2,
                                tap_module_min=0.5,
                                shift_angle=0.0,
                                theta_max=6.28,
                                theta_min=-6.28,
                                active=True,)

        data.append(elm)

    return data


def parse_branches(data_lst: List[List], buses_dict: Dict[int, dev.Bus]):
    data = list()
    for raw in data_lst:
        name = '{0}_{1}_{2}'.format(raw[1], raw[4], raw[6])
        code = '{0}_{1}_{2}'.format(raw[0], raw[3], raw[6])
        bus_f = buses_dict[raw[0]]
        bus_t = buses_dict[raw[3]]
        rate = raw[13]
        r = raw[10]
        x = raw[11]
        b = raw[12]
        elm = dev.Line(bus_from=bus_f,
                       bus_to=bus_t,
                       name=name,
                       idtag=None,
                       code=code,
                       r=r,
                       x=x,
                       b=b,
                       rate=rate,
                       active=True,)

        data.append(elm)

    return data


def parse_dc_lines(data_lst: List[List], buses_dict: Dict[int, dev.Bus], Sbase=100):
    data = list()
    for raw in data_lst:
        name = '{0}_{1}_{2}'.format(raw[1], raw[4], raw[6])
        code = '{0}_{1}_{2}'.format(raw[0], raw[3], raw[6])
        bus_f = buses_dict[raw[0]]
        bus_t = buses_dict[raw[3]]
        rate = raw[14]
        zbase = bus_f.Vnom * bus_f.Vnom / Sbase
        r = raw[11] / zbase
        elm = dev.DcLine(bus_from=bus_f,
                         bus_to=bus_t,
                         name=name,
                         idtag=None,
                         code=code,
                         r=r,
                         rate=rate,
                         active=True,)

        data.append(elm)

    return data


def parse_dc_converters(data_lst: List[List], buses_dict: Dict[int, dev.Bus], dc_buses_dict: Dict[int, dev.Bus]):
    data = list()
    for raw in data_lst:
        name = '{0}_{1}_{2}'.format(raw[1], raw[4], raw[6])
        code = '{0}_{1}_{2}'.format(raw[0], raw[3], raw[6])
        bus_t = buses_dict[raw[0]]
        bus_f = dc_buses_dict[raw[3]]
        rate = 100.0

        elm = dev.VSC(bus_from=bus_f,
                      bus_to=bus_t,
                      name=name,
                      idtag=None,
                      code=code,
                      rate=rate,
                      active=True,)

        data.append(elm)

    return data


def parse_loads(data_lst: List[List], buses_dict: Dict[int, dev.Bus]):
    data = list()
    for raw in data_lst:
        name = '{0}_{1}'.format(raw[1], raw[3])
        code = '{0}_{1}'.format(raw[0], raw[3])
        bus_f = buses_dict[raw[0]]
        P = raw[6]
        Q = raw[7]

        elm = dev.Load(name=name,
                       idtag=None,
                       code=code,
                       P=P,
                       Q=Q,
                       active=True,)
        bus_f.add_device(elm)

        data.append(elm)

    return data


def parse_generators(data_lst: List[List], buses_dict: Dict[int, dev.Bus], bus_volt: Dict[int, float]):
    data = list()
    for raw in data_lst:
        name = '{0}_{1}'.format(raw[1], raw[3])
        code = '{0}_{1}'.format(raw[0], raw[3])
        bus_f = buses_dict[raw[0]]
        Vset = bus_volt[raw[0]]
        P = raw[13]
        Pmax = raw[14]
        Pmin = raw[15]
        Q = raw[16]
        Qmax = raw[17]
        Qmin = raw[18]
        Sbase = raw[19]

        elm = dev.Generator(name=name,
                            idtag=None,
                            code=code,
                            active_power=P,
                            voltage_module=Vset,
                            p_min=Pmin,
                            p_max=Pmax,
                            Qmin=Qmin,
                            Qmax=Qmax,
                            Snom=np.sqrt(Pmax*Pmax+Qmax*Qmax),
                            Sbase=Sbase,
                            active=True, )
        bus_f.add_device(elm)

        data.append(elm)

    return data


class PowerWorldParser:

    def __init__(self, file_name):
        """
        Parse PowerWorld EPC file
        Args:
            file_name: file name or path
        """
        self.parsers = dict()
        self.versions = []

        self.logger = Logger()

        self.file_name = file_name

        self.circuit, self.logger = self.parse_case()

        self.circuit.comments = 'Converted from the PowerWorld .epc file ' \
                                + os.path.basename(file_name) + '\n\n' + str(self.logger)

    def read_and_split(self) -> (List[AnyStr], Dict[AnyStr, AnyStr]):
        """
        Read the text file and split it into sections
        :return: list of sections, dictionary of sections by type
        """

        # make a guess of the file encoding
        detection = chardet.detect(open(self.file_name, "rb").read())

        # open the text file into a variable
        txt = ''
        with open(self.file_name, 'r', encoding=detection['encoding']) as my_file:
            for line in my_file:
                if line[0] != '@':
                    txt += line

        # fix stupid line partition
        txt = txt.replace('/\n', '')

        expected_sections = ['title',
                             'comments',
                             'solution parameters',
                             'substation data',
                             'bus data',
                             'branch data',
                             'transformer data',
                             'generator data',
                             'load data',
                             'shunt data',
                             'svd data',
                             'area data',
                             'zone data',
                             'interface data',
                             'interface branch data',
                             'dc bus data',
                             'dc line data',
                             'dc converter data',
                             'z table data',
                             'gcd data',
                             'transaction data',
                             'owner data',
                             'qtable data',
                             'ba data',
                             'end']

        # find which of the expected are actually there
        present_sections = list()
        for a in expected_sections:
            if a in txt:
                present_sections.append(a)

        # split the text file into sections
        sections_dict = dict()
        for i in range(len(present_sections)-1):
            a = present_sections[i]
            b = present_sections[i + 1]
            if a in txt and b in txt:
                raw_txt = find_between(txt, a, b)
                lines = raw_txt.split('\n')

                if len(lines) > 0:
                    if '[' in lines[0]:
                        new_lines = list()
                        header = lines[0].split(']')[1].split()
                        for j in range(1, len(lines)):
                            line_data = split_line(lines[j])
                            if len(line_data) > 0:
                                new_lines.append(line_data)

                        sections_dict[a] = {'header': header, 'data': new_lines}
                else:
                    sections_dict[a] = {'header': '', 'data': lines}
            else:
                sections_dict[a] = {'header': '', 'data': list()}

        return sections_dict

    def parse_case(self) -> (MultiCircuit, List[AnyStr]):
        """
        EPC power world case

        Returns: MultiCircuit, List[str]
        """
        grid = MultiCircuit()
        logger = Logger()

        data_dict = self.read_and_split()

        if 'substation data' in data_dict.keys():
            substations_dict = parse_substations(data_dict['substation data']['data'])
            grid.substations = list(substations_dict.values())
        else:
            substations_dict = None

        buses_dict, bus_volt = parse_buses(data_dict['bus data']['data'], substations_dict)

        # create devices
        grid.buses = list(buses_dict.values())

        if 'branch data' in data_dict.keys():
            grid.lines = parse_branches(data_dict['branch data']['data'], buses_dict)

        if 'transformer data' in data_dict.keys():
            grid.transformers2w = parse_transformers(data_dict['transformer data']['data'], buses_dict)

        if 'load data' in data_dict.keys():
            parse_loads(data_dict['load data']['data'], buses_dict)

        if 'generator data' in data_dict.keys():
            parse_generators(data_dict['generator data']['data'], buses_dict, bus_volt)

        if 'dc bus data' in data_dict.keys():
            # augments buses_dict and bus_volt
            dc_buses_dict, dc_bus_volt = parse_dc_buses(data_dict['dc bus data']['data'])
            grid.buses += list(dc_buses_dict.values())

            if 'dc line data' in data_dict.keys():
                grid.dc_lines = parse_dc_lines(data_dict['dc line data']['data'], dc_buses_dict)

            if 'dc converter data' in data_dict.keys():
                grid.vsc_devices = parse_dc_converters(data_dict['dc converter data']['data'], buses_dict, dc_buses_dict)

        grid.fill_xy_from_lat_lon()

        return grid, logger


if __name__ == '__main__':

    # f = '/home/santi/Descargas/ACTIVSg500/ACTIVSg500.EPC'
    f = r'C:\Users\SPV86\Downloads\IEEE300\IEEE300Bus.epc'

    parser = PowerWorldParser(f)
    parser.parse_case()

    print()