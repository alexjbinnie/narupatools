# This file is part of narupatools (https://github.com/alexjbinnie/narupatools).
# Copyright (c) Alex Jamieson-Binnie. All rights reserved.
#
# narupatools is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# narupatools is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with narupatools.  If not, see <http://www.gnu.org/licenses/>.
#
# Originally part of the narupa-openmm package.
# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Modified under the terms of the GPL.

"""More performant serialization and deserialization of OpenMM simulations."""

from io import StringIO

from lxml import etree
from lxml.etree import ElementBase
from simtk.openmm import Integrator, System, XmlSerializer
from simtk.openmm.app import PDBFile, Simulation


def _get_single_node(document: ElementBase, tag_name: str) -> ElementBase:
    nodes = document.findall(tag_name)
    if len(nodes) == 0:
        raise IOError(f"No {tag_name} tag defined in the XML.")
    if len(nodes) != 1:
        raise IOError(f"More than one {tag_name} tag defined in the XML.")
    return nodes[0]


def serialize_simulation(simulation: Simulation) -> str:
    """
    Serialize an OpenMM simulation to an XML file.

    This is more performant than the equivalent method provided by Narupa, as it uses a
    faster XML library.

    :param simulation: OpenMM simulation to serialize.
    :return: String that contains XML serialized simulation.
    """
    positions = simulation.context.getState(getPositions=True).getPositions()
    pdb_content = StringIO()
    PDBFile.writeFile(simulation.topology, positions, pdb_content)
    pdb_content = pdb_content.getvalue()

    # Extract the system
    system_xml_str = XmlSerializer.serialize(simulation.system)
    system_document = etree.fromstring(system_xml_str)

    # Extract the integrator
    integrator_xml_str = XmlSerializer.serialize(simulation.integrator)
    integrator_document = etree.fromstring(integrator_xml_str)

    # Combine the element in a single
    root = etree.Element("OpenMMSimulation")
    pdb_node = etree.SubElement(root, "pdb")
    pdb_node.text = pdb_content
    root.append(system_document)
    root.append(integrator_document)

    return etree.tostring(root, encoding="unicode", pretty_print=True)


def deserialize_simulation(contents: str) -> Simulation:
    """
    Deserialize an XML string to an OpenMM simulation.

    This is more performant than the equivalent method provided by Narupa, as it uses a
    faster XML library.

    :param contents: Contents of the XML file.
    :return: OpenMM simulation read from the provided string.
    """
    document = etree.fromstring(contents)

    pdb_node = _get_single_node(document, "pdb")
    pdb_io = StringIO(pdb_node.text)
    pdb = PDBFile(pdb_io)

    system_node = _get_single_node(document, "System")
    system_content = etree.tostring(system_node, encoding="unicode", method="xml")
    system: System = XmlSerializer.deserialize(system_content)  # type: ignore

    integrator_node = _get_single_node(document, "Integrator")
    integrator_content = etree.tostring(
        integrator_node, encoding="unicode", method="xml"
    )
    integrator: Integrator = XmlSerializer.deserialize(
        integrator_content
    )  # type: ignore

    simulation = Simulation(
        topology=pdb.topology,
        system=system,
        integrator=integrator,
    )
    simulation.context.setPositions(pdb.positions)
    return simulation
