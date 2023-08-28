from typing import List
from enum import Enum
import datetime
import pandas as pd


class DataLogSeverity(Enum):
    Error = 'Error'
    Warning = 'Warning'
    Information = 'Information'
    Divergence = 'Divergence'

    def __str__(self):
        return self.value

    def __repr__(self):
        return str(self)

    @staticmethod
    def argparse(s):
        try:
            return DataLogSeverity[s]
        except KeyError:
            return s


class DataLogEntry:

    def __init__(self, msg="", severity: DataLogSeverity = DataLogSeverity.Information, device="", device_class="",
                 property_name='', value="", expected_value="", comment=""):
        self.time = "{date:%H:%M:%S}".format(date=datetime.datetime.now())  # might use %Y/%m/%d %H:%M:%S
        self.msg: str = str(msg)
        self.severity: DataLogSeverity = severity
        self.device = device
        self.device_class = device_class
        self.property_name = property_name
        self.value = value
        self.expected_value = str(expected_value)
        self.comment = comment

    def to_list(self):
        return [self.time, self.severity.value, self.msg, self.device, self.device_class,
                self.property_name, self.value, self.expected_value, self.comment]

    def __str__(self):
        return "{0} {1}: {2} {3} {4} {5} {6} {7} {8}".format(self.time,
                                                             self.severity.value,
                                                             self.msg,
                                                             self.device,
                                                             self.device_class,
                                                             self.property_name,
                                                             self.value,
                                                             self.expected_value,
                                                             self.comment)


class DataLogger:

    def __init__(self):

        self.entries: List[DataLogEntry] = list()

        self.debug_entries: List[str] = list()

    def get_message(self):
        """
        Get a diagnostic message
        :return: message
        """
        n_warning = 0
        n_error = 0
        n_info = 0

        for m in self.entries:
            if m.severity == DataLogSeverity.Error:
                n_error += 1
            elif m.severity == DataLogSeverity.Warning:
                n_warning += 1
            elif m.severity == DataLogSeverity.Information:
                n_info += 1

        return "There were {} errors, {} warnings and {} info logs.".format(n_error, n_warning, n_info)

    def add_debug(self, *args):
        self.debug_entries.append(" ".join([str(x) for x in args]))

    def append(self, txt: str):
        """

        :param txt:
        :return:
        """
        self.entries.append(DataLogEntry(txt))

    def has_logs(self):
        return len(self.entries) > 0

    def add_info(self, msg: str, device="", device_class="", device_property='', value="", expected_value=""):
        """

        :param msg:
        :param device:
        :param value:
        :param expected_value
        :return:
        """
        self.entries.append(DataLogEntry(msg, DataLogSeverity.Information, device, device_class, device_property,
                                         str(value), str(expected_value)))

    def add_warning(self, msg, device="", device_class="", device_property='', value="", expected_value=""):
        """

        :param msg:
        :param device:
        :param device_class:
        :param property:
        :param value:
        :param expected_value:
        :return:
        """
        self.entries.append(
            DataLogEntry(msg, DataLogSeverity.Warning, device, device_class, device_property, str(value),
                         str(expected_value)))

    def add_error(self, msg, device="", device_class="", device_property='', value="", expected_value="", comment=""):
        """

        :param msg:
        :param device:
        :param device_class:
        :param device_property:
        :param value:
        :param expected_value:
        :param comment:
        :return:
        """
        self.entries.append(DataLogEntry(msg=msg,
                                         severity=DataLogSeverity.Error,
                                         device=device,
                                         device_class=device_class,
                                         property_name=device_property,
                                         value=str(value),
                                         expected_value=str(expected_value),
                                         comment=comment))

    def add_divergence(self, msg, device="", device_class="", device_property='', value=0, expected_value=0, tol=1e-6):
        """

        :param msg:
        :param device:
        :param device_class:
        :param device_property:
        :param value:
        :param expected_value:
        :param tol:
        :return:
        """

        if abs(value - expected_value) > tol:
            self.entries.append(DataLogEntry(msg=msg,
                                             severity=DataLogSeverity.Divergence,
                                             device=device,
                                             device_class=device_class,
                                             property_name=device_property,
                                             value=str(value),
                                             expected_value=str(expected_value)))

    def add(self, msg, severity: DataLogSeverity = DataLogSeverity.Error, device="", device_class="",
            device_property='', value="", expected_value=""):
        """

        :param msg:
        :param severity:
        :param device:
        :param device_class:
        :param device_property:
        :param value:
        :param expected_value:
        :return:
        """
        self.entries.append(DataLogEntry(msg=msg,
                                         severity=severity,
                                         device=device,
                                         device_class=device_class,
                                         property_name=device_property,
                                         value=str(value),
                                         expected_value=str(expected_value)))

    def to_dict(self):
        """
        Get the logs sorted by severity and message
        :return: Dictionary[Dictionary[List[time, device, value, expected value]]]
        """
        by_severity = dict()

        for e in self.entries:

            if e.severity.value not in by_severity.keys():
                by_severity[e.severity.value] = dict()

            by_msg = by_severity[e.severity.value]

            if e.msg in by_msg.keys():
                # add msg to existing msg list
                by_msg[e.msg].append((e.time, e.device, e.device_class, e.property_name,
                                      e.value, e.expected_value, e.comment))
            else:
                # add msg entry for the first time
                by_msg[e.msg] = [(e.time, e.device, e.device_class, e.property_name,
                                  e.value, e.expected_value, e.comment)]

        return by_severity

    def to_df(self):
        """
        Get DataFrame
        :return:
        """
        data = [e.to_list() for e in self.entries]
        df = pd.DataFrame(data=data, columns=['Time', 'Severity', 'Message', 'Device', 'Property',
                                              'Class', 'Value', 'Expected value', 'Comment'])
        df.set_index('Time', inplace=True)
        return df

    def to_csv(self, fname):
        self.to_df().to_csv(fname)

    def to_xlsx(self, fname, sheet_name='Logs'):
        self.to_df().to_excel(fname, sheet_name=sheet_name)

    def __str__(self):

        val = ''
        for e in self.entries:
            val += str(e) + '\n'
        return val

    def __getitem__(self, key):
        """
        get [index] implementation
        :param key: integer
        :return: message, severity
        """
        return self.entries[key]

    def __setitem__(self, idx, value):
        """
        set [index] implementation
        :param idx: integer
        :param value: string message
        :return: Nothing
        """
        self.entries[idx] = value

    def __iadd__(self, other: "Logger"):
        """
        += implementation
        :param other:
        :return:
        """

        if other is not None:
            self.entries += other.entries
        return self

    def __len__(self):
        return len(self.entries)

    def size(self):
        return len(self.entries)
