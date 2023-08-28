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
from GridCal.Engine.IO.cim.cgmes_2_4_15.devices.identified_object import IdentifiedObject, cgmesProfile


class FullModel(IdentifiedObject):

    def __init__(self, rdfid, tpe):
        IdentifiedObject.__init__(self, rdfid, tpe)

        self.scenarioTime = ''
        self.created = ''
        self.version = ''
        self.profile = ''
        self.modelingAuthoritySet = ''
        self.DependentOn = ''
        self.longDependentOnPF = ''
        self.Supersedes = ''

        self.register_property(
            name='scenarioTime',
            class_type=str,
            description="scenarioTime.",
            profiles=[cgmesProfile.TP_BD, cgmesProfile.DL, cgmesProfile.SSH, cgmesProfile.EQ,
                      cgmesProfile.DY, cgmesProfile.TP, cgmesProfile.EQ_BD, cgmesProfile.GL,
                      cgmesProfile.SV])

        self.register_property(
            name='created',
            class_type=str,
            description="Creation date.",
            profiles=[cgmesProfile.TP_BD, cgmesProfile.DL, cgmesProfile.SSH, cgmesProfile.EQ,
                      cgmesProfile.DY, cgmesProfile.TP, cgmesProfile.EQ_BD, cgmesProfile.GL,
                      cgmesProfile.SV])

        self.register_property(
            name='version',
            class_type=int,
            description="version.",
            profiles=[cgmesProfile.TP_BD, cgmesProfile.DL, cgmesProfile.SSH, cgmesProfile.EQ,
                      cgmesProfile.DY, cgmesProfile.TP, cgmesProfile.EQ_BD, cgmesProfile.GL,
                      cgmesProfile.SV])

        self.register_property(
            name='profile',
            class_type=str,
            description="profile.",
            profiles=[cgmesProfile.TP_BD, cgmesProfile.DL, cgmesProfile.SSH, cgmesProfile.EQ,
                      cgmesProfile.DY, cgmesProfile.TP, cgmesProfile.EQ_BD, cgmesProfile.GL,
                      cgmesProfile.SV])

        self.register_property(
            name='modelingAuthoritySet',
            class_type=str,
            description="modelingAuthoritySet",
            profiles=[cgmesProfile.TP_BD, cgmesProfile.DL, cgmesProfile.SSH, cgmesProfile.EQ,
                      cgmesProfile.DY, cgmesProfile.TP, cgmesProfile.EQ_BD, cgmesProfile.GL,
                      cgmesProfile.SV])

        self.register_property(
            name='DependentOn',
            class_type=str,
            description="DependentOn.",
            profiles=[cgmesProfile.TP_BD, cgmesProfile.DL, cgmesProfile.SSH, cgmesProfile.EQ,
                      cgmesProfile.DY, cgmesProfile.TP, cgmesProfile.EQ_BD, cgmesProfile.GL,
                      cgmesProfile.SV])

        self.register_property(
            name='longDependentOnPF',
            class_type=str,
            description="longDependentOnPF.",
            profiles=[cgmesProfile.TP_BD, cgmesProfile.DL, cgmesProfile.SSH, cgmesProfile.EQ,
                      cgmesProfile.DY, cgmesProfile.TP, cgmesProfile.EQ_BD, cgmesProfile.GL,
                      cgmesProfile.SV])

        self.register_property(
            name='Supersedes',
            class_type=str,
            description="Supersedes.",
            profiles=[cgmesProfile.TP_BD, cgmesProfile.DL, cgmesProfile.SSH, cgmesProfile.EQ,
                      cgmesProfile.DY, cgmesProfile.TP, cgmesProfile.EQ_BD, cgmesProfile.GL,
                      cgmesProfile.SV])
