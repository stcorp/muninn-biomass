import os
import re
import zipfile
from datetime import datetime
import xml.etree.ElementTree as ET
from typing import Optional, Callable

from muninn.geometry import Point, LinearRing, Polygon, MultiPolygon
from muninn.schema import Mapping, Text, Integer, Boolean, Timestamp
from muninn.util import copy_path
from muninn import Struct
from muninn import Error as muninnError


# Namespaces

class BIOMASSNamespace(Mapping):
    mission_phase = Text(index=True, optional=True)  # C (Commissioning), I (Interferometric), T (Tomographic)
    global_coverage = Integer(index=True, optional=True)  # 1 (TOM), 1-6 (INT)
    major_cycle = Integer(index=True, optional=True)  # 1-7 (TOM/INT)
    repeat_cycle = Integer(index=True, optional=True)  # 1-8 or 0 for drift (not used for L2 products)
    track_number = Integer(index=True, optional=True)  # relative orbit within repeat cycle
    frame_number = Integer(index=True, optional=True)  # only applicable for L1/L2A
    slice_number = Integer(index=True, optional=True)  # only applicable for L0
    tile = Text(index=True, optional=True)  # [N|S]aa[E|W]bbb, only applicable for L2/L3, space-separated list for L2A
    basin = Text(index=True, optional=True)  # Basin identifier, only applicable for L2/L3, space-separated list for L2A
    baseline = Integer(index=True, optional=True)  # seperate sequence for each L0/L1/L2 x COM/TOM/INT combination
    compact_creation_date = Text(index=True, optional=True)  # base 36 ([0-9A-Z]) of [s] since 2000-01-01T00:00:00
    orbit_number = Integer(index=True, optional=True)
    orbit_drift_flag = Boolean(index=True, optional=True)  # drift (datatake inside Spacecraft Repositioning Manoeuvre)
    orbit_direction = Text(index=True, optional=True)  # ascending / descending
    swath = Text(index=True, optional=True)  # S1, S2, S3, RO, EC, AC
    type = Text(index=True, optional=True)  # RAW, SCS, DGM, STA, GN, FD, FH, AGB, ORB, ATT, INS, PP0, PP1, PP2, PP3,
                                            # PPS, TEC, GMF, TRP
    processing_level = Text(index=True, optional=True)  # RAW, AUX, L0, L1A, L1B, L1C, L2A, L2B, L3
    datatake_id = Text(index=True, optional=True)  # decimal number - can be space-separated list for RAWS products
    instrument_conf_id = Integer(index=True, optional=True)
    processing_center = Text(index=True, optional=True)
    processor_name = Text(index=True, optional=True)
    processor_version = Text(index=True, optional=True)
    acquisition_station = Text(index=True, optional=True)
    acquisition_date = Timestamp(index=True, optional=True)


def namespaces():
    return ["biomass"]


def namespace(namespace_name):
    return BIOMASSNamespace


# Product types

L0_PRODUCT_TYPES = [
    "S1_RAW__0S",
    "S2_RAW__0S",
    "S3_RAW__0S",
    "S1_RAW__0M",
    "S2_RAW__0M",
    "S3_RAW__0M",
    "RO_RAW__0S",
    "EC_RAW__0S",
    "AC_RAW__0A",
]
L1_PRODUCT_TYPES = [
    "S1_SCS__1S",
    "S2_SCS__1S",
    "S3_SCS__1S",
    "S1_DGM__1S",
    "S2_DGM__1S",
    "S3_DGM__1S",
    "S1_SCS__1M",
    "S2_SCS__1M",
    "S3_SCS__1M",
    "S1_STA__1S",
    "S2_STA__1S",
    "S3_STA__1S",
    "S1_STA__1M",
    "S2_STA__1M",
    "S3_STA__1M",
    "RO_SCS__1S",
]
L2A_PRODUCT_TYPES = [
    "FP_GN__L2A",
    "FP_FD__L2A",
    "FP_FH__L2A",
]
L2B_L3_PRODUCT_TYPES = [
    "FP_FD__L2B",
    "FP_FH__L2B",
    "FP_AGB_L2B",
    "FP_FD__L3T",
    "FP_FH__L3T",
    "FP_AGB_L3T",
]
ADF_PRODUCT_TYPES = [
    "AUX_ORB___",
    "AUX_ATT___",
    "AUX_INS___",
    "AUX_PP0___",
    "AUX_PP1___",
    "AUX_PPS___",
    "AUX_PP2_2A",
    "AUX_PP2_FH",
    "AUX_PP2_FD",
    "AUX_PP2_AB",
    "AUX_PP3___",
    "AUX_TEC___",
    "AUX_GMF___",
]
BASELINE_CFG_PRODUCT_TYPES = [
    "CFG_PF_BC_",
]
REPORT_PRODUCT_TYPES = [
    "REP_NOB_IN",
    "REP_NOBOUT",
    "REP_PROD__",
]
QUALITY_DISCLAIMER_PRODUCT_TYPES = [
    "REP_DISCLM",
]
TRANSPONDER_EOF_FILE_TYPES = [
    "CFG_TRPCCF",
    "REP_TRPCRF",
    "REP_TRPSCR",
    "LOG_TRPSLF",
    "LOG_TRPDLF",
]
TRANSPONDER_REPORT_TYPES = [
    "REP_TRPACQ",
    "REP_TRPICR",
    "REP_TRPCAL",
]
MISSION_PLANNING_EOF_TYPES = [
    "MPL_TRPINP",
    "MPL_XBSPLN",
]
MISSION_PLANNING_EOF_WITH_BASELINE_TYPES = [
    "MPL_PPF_1_",
    "MPL_PPF_2_",
    "MPL_PPF_3_",
    "MPL_PPF_4_",
    "MPL_PPF_5_",
    "MPL_PPF_6_",
    "MPL_PPF_7_",
    "MPL_PPF_8_",
    "MPL_ORBREF",
    "MPL_ORBPRE",
]
MISSION_PLANNING_EOF_WITH_ACQ_TYPES = [
    "MPL_XBSACQ",
]
MISSION_PLANNING_EOF_WITH_PHASE_TYPES = [
    "MPL_MISSTL",
]


def compress(paths, target_filepath, compresslevel=None):
    if compresslevel is None:
        compression = zipfile.ZIP_STORED
    else:
        compression = zipfile.ZIP_DEFLATED
    with zipfile.ZipFile(target_filepath, "x", compression, compresslevel=compresslevel, strict_timestamps=False) \
            as archive:
        for path in paths:
            rootlen = len(os.path.dirname(path)) + 1
            if os.path.isdir(path):
                for base, dirs, files in os.walk(path):
                    for file in files:
                        fn = os.path.join(base, file)
                        archive.write(fn, fn[rootlen:])
            else:
                archive.write(path, path[rootlen:])


def parse_boolean(val):
    val = val.lower()
    if val == "true":
        res = True
    elif val == "false":
        res = False
    else:
        raise ValueError(f"Invalid value to convert to bool: {val}")
    return res


class BiomassBaseProduct(object):
    # filename_base_pattern is the pattern for the filename excluding any extension (and without trailing $)
    # extension can be "" (for a directory), None (for multifile products), or set to a specific extension (e.g. ".EOF")
    def __init__(self, product_type: str, filename_base_pattern: str = None, extension: str = None,
                 zipped: Optional[bool] = None):
        self.is_multi_file_product = False
        self.product_type = product_type
        self.extension = extension
        self.filename_pattern = filename_base_pattern
        if zipped is None:  # "None" means flexible zip handling
            if extension == "":
                self.filename_pattern += r"(\.ZIP)?$"
            elif extension is None:
                self.filename_pattern += r"(\.ZIP$)?"
            else:
                self.filename_pattern += r"(%s|\.ZIP)$" % re.escape(extension)
        elif zipped:
            self.filename_pattern += r"\.ZIP$"
        else:
            if self.extension is not None:
                self.filename_pattern += r"%s$" % re.escape(extension)

    @property
    def use_enclosing_directory(self):
        return False

    @property
    def namespaces(self):
        return ["biomass"]

    @staticmethod
    def archive_path(attributes):
        return os.path.join(
            attributes.core.product_type,
            attributes.core.validity_start.strftime("%Y"),
            attributes.core.validity_start.strftime("%m"),
            attributes.core.validity_start.strftime("%d"),
        )

    def is_zipped(self, filepath):
        return filepath.endswith(".ZIP")

    def parse_filename(self, filename):
        match = re.match(self.filename_pattern, os.path.basename(filename))
        if match:
            return match.groupdict()
        return None

    def identify(self, paths):
        if len(paths) != 1:
            return False
        return re.match(self.filename_pattern, os.path.basename(paths[0])) is not None

    def analyze(self, paths, filename_only=False):
        file_path = paths[0]
        file_name = os.path.basename(file_path)
        file_name_attrs = self.parse_filename(file_name)

        properties = Struct()

        core = properties.core = Struct()
        core.product_name = os.path.splitext(file_name)[0]
        properties.biomass = Struct()

        self._set_properties_from_filename(file_name_attrs, properties)

        return properties

    def _set_properties_from_filename(self, file_name_attrs, properties):
        core = properties.core
        core.validity_start = datetime.strptime(file_name_attrs["validity_start"], "%Y%m%dT%H%M%S")
        if file_name_attrs["validity_stop"] == "99999999T999999":
            core.validity_stop = datetime.max
        else:
            core.validity_stop = datetime.strptime(file_name_attrs["validity_stop"], "%Y%m%dT%H%M%S")

    def read_xml_component(self, filepath, component_path=None):
        # filepath: Path given as input to the analyze() function
        # component_path: Path of the specific component to be read.

        # Open XML file (zipped or not) and return root element
        if self.is_zipped(filepath):
            if component_path is None:
                component_path = os.path.splitext(os.path.basename(filepath))[0]
                if self.extension is not None:
                    component_path += self.extension
            else:
                if not self.is_multi_file_product:
                    component_path = os.path.join(os.path.splitext(os.path.basename(filepath))[0], component_path)
            with zipfile.ZipFile(filepath) as zproduct:
                with zproduct.open(component_path) as file:
                    return ET.parse(file).getroot()
        else:
            if component_path is not None:
                filepath = os.path.join(filepath, component_path)
            with open(filepath) as file:
                return ET.parse(file).getroot()

    @staticmethod
    def _set_property(properties: Struct, property_name: str, root: ET.Element, path: str, ns: dict, parse: Callable,
                      mandatory: Optional[bool] = False) -> None:
        """
        Helper method to handle assigning properties

        Args
            properties (Struct): Muninn namespace element where the property must be set
            property_name (str): Name of the property to set (as given in the namespace)
            root (ET.Element): Root of the ET.ElementTree
            path (str): Absolute path of the element to retrieve the property from
            ns (dict): Contains namespace information necessary for ET.iterfind
            parse (function): Function to convert the string read from the ET element to the desired type.
                                If the function requires additional parameters, provide a lambda function
                                that includes the parameters.
            mandatory (bool): Flag that controls if an error should be thrown if this property cannot be extracted:
                                False -> No error. Defaults to False.
        """
        property_iter = root.iterfind(path, namespaces=ns)
        list_values = []
        for tmp in property_iter:
            if tmp is not None:
                list_values.append(tmp.text)

        property_value = " ".join(list_values)
        if property_value != "":
            try:
                value = parse(property_value)
                if value is not None:
                    # Value is only None if the metadata was signified as not available (e.g. "NA" or "___").
                    # In that case, the metadata is ignored. If not, it is stored into properties.
                    properties[property_name] = value
            except ValueError as e:
                raise muninnError(f"could not extract metadata {property_name} ({e})")
            except IndexError:
                raise muninnError(f"could not extract metadata {property_name} (incorrect value: {property_value})")

        elif mandatory and not hasattr(properties, property_name):
            raise muninnError(f"could not extract mandatory metadata {property_name} (metadata missing)")

    def _analyze_mph(self, root, properties):
        ns = {
            "bio": "http://earth.esa.int/biomass/1.0",
            "eop": "http://www.opengis.net/eop/2.1",
            "gml": "http://www.opengis.net/gml/3.2",
            "om": "http://www.opengis.net/om/2.0",
            "ows": "http://www.opengis.net/ows/2.0",
            "sar": "http://www.opengis.net/sar/2.1",
            "xlink": "http://www.w3.org/1999/xlink",
            "xsi": "http://www.w3.org/2001/XMLSchema-instance",
        }

        biomass = properties.biomass
        core = properties.core
        product_type = self.product_type

        mandatory = {
            "swath": product_type in L0_PRODUCT_TYPES + L1_PRODUCT_TYPES + L2A_PRODUCT_TYPES,
            "orbit_number": product_type in L0_PRODUCT_TYPES + L1_PRODUCT_TYPES,
            "orbit_drift_flag": (product_type in L0_PRODUCT_TYPES + L1_PRODUCT_TYPES) and
                                (product_type != "AC_RAW__0A"),
            "orbit_direction": product_type in L0_PRODUCT_TYPES + L1_PRODUCT_TYPES+L2A_PRODUCT_TYPES,
            "datatake_id": product_type in L0_PRODUCT_TYPES + L1_PRODUCT_TYPES + ["AUX_ATT___", "AUX_ORB___"],
            "instrument_conf_id": (product_type in L0_PRODUCT_TYPES + L1_PRODUCT_TYPES) and
                                  (product_type != "AC_RAW__0A"),
            "track_number": product_type in L0_PRODUCT_TYPES + L1_PRODUCT_TYPES + L2A_PRODUCT_TYPES,
            "slice_number": product_type in L0_PRODUCT_TYPES,
            "frame_number": product_type in L1_PRODUCT_TYPES + L2A_PRODUCT_TYPES,
            "mission_phase": product_type in (L0_PRODUCT_TYPES + L1_PRODUCT_TYPES + L2A_PRODUCT_TYPES +
                                              L2B_L3_PRODUCT_TYPES),
            "global_coverage": product_type in (L0_PRODUCT_TYPES + L1_PRODUCT_TYPES + L2A_PRODUCT_TYPES +
                                                L2B_L3_PRODUCT_TYPES),
            "major_cycle": product_type in L0_PRODUCT_TYPES + L1_PRODUCT_TYPES + L2A_PRODUCT_TYPES,
            "repeat_cycle": product_type in L0_PRODUCT_TYPES + L1_PRODUCT_TYPES,
            "tile": product_type in L2A_PRODUCT_TYPES + L2B_L3_PRODUCT_TYPES,
            "basin": product_type in L2B_L3_PRODUCT_TYPES,
            "acquisition_station": False,  # Mandatory for RAW products, but these are not yet implemented
            "acquisition_date": False,  # Mandatory for RAW products, but these are not yet implemented
        }

        pathSensor = "./om:procedure/eop:EarthObservationEquipment/eop:sensor/eop:Sensor/"
        property_name = "swath"
        self._set_property(biomass, property_name, root, pathSensor + "eop:swathIdentifier", ns, str,
                           mandatory[property_name])

        pathAcquisition = "./om:procedure/eop:EarthObservationEquipment/eop:acquisitionParameters/bio:Acquisition/"
        property_name = "orbit_number"
        self._set_property(biomass, property_name, root, pathAcquisition + "eop:orbitNumber", ns, int,
                           mandatory[property_name])
        property_name = "orbit_drift_flag"
        self._set_property(biomass, property_name, root, pathAcquisition + "bio:orbitDriftFlag", ns, parse_boolean,
                           mandatory[property_name])
        property_name = "orbit_direction"
        self._set_property(biomass, property_name, root, pathAcquisition + "eop:orbitDirection", ns, str,
                           mandatory[property_name])
        property_name = "datatake_id"
        self._set_property(biomass, property_name, root, pathAcquisition + "bio:dataTakeID", ns, str,
                           mandatory[property_name])
        property_name = "instrument_conf_id"
        self._set_property(biomass, property_name, root, pathAcquisition + "bio:instrumentConfID", ns, int,
                           mandatory[property_name])
        property_name = "track_number"
        self._set_property(biomass, property_name, root, pathAcquisition + "eop:wrsLongitudeGrid", ns, int,
                           mandatory[property_name])
        if self.product_type in L0_PRODUCT_TYPES:
            # For L0 product, this contains the slice number
            property_name = "slice_number"
            self._set_property(biomass, property_name, root, pathAcquisition + "eop:wrsLatitudeGrid", ns,
                               lambda x: None if "_" in x else int(x),
                               mandatory[property_name])  # Not available value ("___") is ignored
        elif self.product_type in L1_PRODUCT_TYPES+L2A_PRODUCT_TYPES:
            # For L1 and L2A, this contains the frame number
            property_name = "frame_number"
            self._set_property(biomass, property_name, root, pathAcquisition + "eop:wrsLatitudeGrid", ns, int,
                               mandatory[property_name])
        property_name = "mission_phase"
        self._set_property(biomass, property_name, root, pathAcquisition + "bio:missionPhase", ns, lambda x: x[0],
                           mandatory[property_name])  # Only the first letter is taken
        property_name = "global_coverage"
        self._set_property(biomass, property_name, root, pathAcquisition + "bio:globalCoverageID", ns,
                           lambda x: None if x == "NA" else int(x), mandatory[property_name])  # NA is ignored
        property_name = "major_cycle"
        self._set_property(biomass, property_name, root, pathAcquisition + "bio:majorCycleID", ns,
                           lambda x: None if x == "NA" else int(x), mandatory[property_name])  # NA is ignored
        property_name = "repeat_cycle"
        self._set_property(biomass, property_name, root, pathAcquisition + "bio:repeatCycleID", ns,
                           lambda x: 0 if x == "DR" else (None if x == "NA" else int(x)),
                           mandatory[property_name])  # Drift is represented by a 0 value. NA is ignored
        property_name = "tile"
        self._set_property(biomass, property_name, root, pathAcquisition + "bio:tileID", ns, str,
                           mandatory[property_name])
        property_name = "basin"
        self._set_property(biomass, property_name, root, pathAcquisition + "bio:basinID", ns, str,
                           mandatory[property_name])

        pathProcInfo = "./eop:metaDataProperty/bio:EarthObservationMetaData/eop:processing/bio:ProcessingInformation/"
        self._set_property(biomass, "processing_level", root, pathProcInfo + "eop:processingLevel", ns,
                           lambda x: x.split(" ")[1], True)  # Extract LEVEL from "other: LEVEL" in MPH
        self._set_property(biomass, "processing_center", root, pathProcInfo + "eop:processingCenter", ns, str, True)
        self._set_property(biomass, "processor_name", root, pathProcInfo + "eop:processorName", ns, str, True)
        self._set_property(biomass, "processor_version", root, pathProcInfo + "eop:processorVersion", ns, str, True)
        self._set_property(core, "creation_date", root, pathProcInfo + "eop:processingDate", ns,
                           lambda x: datetime.strptime(x, "%Y-%m-%dT%H:%M:%SZ"), True)

        pathDownlinkInfo = "./eop:metaDataProperty/bio:EarthObservationMetaData/eop:downlinkedTo/" \
            "eop:DownlinkInformation/"
        property_name = "acquisition_station"
        self._set_property(biomass, property_name, root, pathDownlinkInfo + "eop:acquisitionStation", ns, str,
                           mandatory[property_name])
        property_name = "acquisition_date"
        self._set_property(biomass, property_name, root, pathDownlinkInfo + "eop:acquisitionDate", ns,
                           lambda x: datetime.strptime(x, "%Y-%m-%dT%H:%M:%SZ"), mandatory[property_name])

        # Extract more accurate versions of validity_start and validity_stop
        pathValidityPeriod = "./om:validTime/gml:TimePeriod/"
        self._set_property(core, "validity_start", root, pathValidityPeriod + "gml:beginPosition", ns,
                           lambda x: datetime.strptime(x, "%Y-%m-%dT%H:%M:%S.%fZ"), True)
        self._set_property(core, "validity_stop", root, pathValidityPeriod + "gml:endPosition", ns,
                           lambda x: datetime.max if x == "9999-99-99T99:99:99.999Z" else
                           datetime.strptime(x, "%Y-%m-%dT%H:%M:%S.%fZ"), True)

        footprint = root.find("./om:featureOfInterest/eop:Footprint", namespaces=ns)
        if footprint is not None:
            polygons = footprint.findall("./eop:multiExtentOf/gml:MultiSurface/gml:surfaceMember/gml:Polygon",
                                         namespaces=ns)
            geometries = []
            for polygon in polygons:
                coord = polygon.find("./gml:exterior/gml:LinearRing/gml:posList", namespaces=ns).text.split(" ")
                geometry = Polygon([LinearRing([Point(float(lon), float(lat)) for lat, lon in
                                    zip(coord[0::2], coord[1::2])])])
                geometries.append(geometry)
            if len(geometries) == 1:
                core.footprint = geometries[0]
            elif len(geometries) > 1:
                core.footprint = MultiPolygon(geometries)

    def _analyze_eof_header(self, root, properties):
        core = properties.core
        biomass = properties.biomass

        # Differentiate between Earth_Observation_Header and Earth_Explorer_Header
        if "Earth_Explorer_File" in root.tag:
            start_node_path = "./Earth_Explorer_Header/"
            # The default namespace changes between some files, so it is extracted manually
            ns = {
                "xsi": "http://www.w3.org/2001/XMLSchema-instance",
                "": root.tag.split("}")[0].split("{")[1] if ("{" in root.tag and "}" in root.tag) else ""
            }
        elif "Earth_Observation_File" in root.tag:
            start_node_path = "./Earth_Observation_Header/"
            ns = {
                "xsi": "http://www.w3.org/2001/XMLSchema-instance",
                "": root.tag.split("}")[0].split("{")[1]  # Same method as above for consistency's sake
            }
        else:
            # Case of a HDR file which starts directly by a "Earth_Observation_Header" element
            start_node_path = "./"
            ns = {
                "xsi": "http://www.w3.org/2001/XMLSchema-instance",
                "": root.tag.split("}")[0].split("{")[1]  # Same method as above for consistency's sake
            }

        # Extract appropriate metadata
        pathSource = start_node_path+"Fixed_Header/Source/"
        self._set_property(biomass, "processing_center", root, pathSource + "System", ns, str, True)
        self._set_property(biomass, "processor_name", root, pathSource + "Creator", ns, str, True)
        self._set_property(biomass, "processor_version", root, pathSource + "Creator_Version", ns, str, True)
        self._set_property(core, "creation_date", root, pathSource + "Creation_Date", ns,
                           lambda x: datetime.strptime(x, "UTC=%Y-%m-%dT%H:%M:%S"), True)

        # Also extract validity_start and validity_stop if not already extracted from the filename
        pathValidityPeriod = start_node_path+"Fixed_Header/Validity_Period/"
        self._set_property(core, "validity_start", root, pathValidityPeriod + "Validity_Start", ns,
                           lambda x: datetime.strptime(x, "UTC=%Y-%m-%dT%H:%M:%S"), True)
        self._set_property(core, "validity_stop", root, pathValidityPeriod + "Validity_Stop", ns,
                           lambda x: datetime.max if x == "UTC=9999-99-99T99:99:99" else
                           datetime.strptime(x, "UTC=%Y-%m-%dT%H:%M:%S"), True)

    def export_zip(self, archive, properties, target_path, paths):
        if self.is_zipped(paths[0]):
            assert len(paths) == 1, "zipped product should be a single file"
            copy_path(paths[0], target_path)
            return os.path.join(target_path, os.path.basename(paths[0]))
        target_filepath = os.path.join(os.path.abspath(target_path), properties.core.physical_name)
        if self.extension:
            target_filepath = target_filepath[:-len(self.extension)]
        target_filepath += ".ZIP"
        compress(paths, target_filepath, compresslevel=1)
        return target_filepath


class FrameBasedDataProduct(BiomassBaseProduct):
    def __init__(self, product_type, zipped=None):
        pattern = [
            r"^BIO",
            product_type,
            r"(?P<validity_start>\d{8}T\d{6})",
            r"(?P<validity_stop>\d{8}T\d{6})",
            r"(?P<mission_phase>[CIT])",
            r"G(?P<global_coverage>[_\d]{2})",
            r"M(?P<major_cycle>[_\d]{2})",
            r"C(?P<repeat_cycle>([_\d]{2}|DR))",
            r"T(?P<track_number>[_\d]{3})",
            r"F(?P<frame_number>[_\d]{3})",
            r"(?P<baseline>\d{2})",
            r"(?P<compact_creation_date>\w{6})",
        ]
        super().__init__(product_type, filename_base_pattern=r"_".join(pattern), extension="", zipped=zipped)

    def _set_properties_from_filename(self, file_name_attrs, properties):
        super()._set_properties_from_filename(file_name_attrs, properties)

        biomass = properties.biomass

        swath = self.product_type.split("_")[0]
        if swath != "FP":  # for L2A we can only derive the swath from the MPH
            biomass.swath = swath
        biomass.type = self.product_type.split("_")[1]
        biomass.mission_phase = file_name_attrs["mission_phase"]
        biomass.global_coverage = int(file_name_attrs["global_coverage"])
        biomass.major_cycle = int(file_name_attrs["major_cycle"])
        if "repeat_cycle" in file_name_attrs and file_name_attrs["repeat_cycle"] != "__":
            if file_name_attrs["repeat_cycle"] == "DR":
                biomass.repeat_cycle = 0
            else:
                biomass.repeat_cycle = int(file_name_attrs["repeat_cycle"])
        if "track_number" in file_name_attrs and file_name_attrs["track_number"] != "___":
            biomass.track_number = int(file_name_attrs["track_number"])
        if "frame_number" in file_name_attrs and file_name_attrs["frame_number"] != "___":
            biomass.frame_number = int(file_name_attrs["frame_number"])
        biomass.baseline = int(file_name_attrs["baseline"])
        biomass.compact_creation_date = file_name_attrs["compact_creation_date"]

    def analyze(self, paths, filename_only=False):
        # Get general BIOMASS product info
        properties = super().analyze(paths, filename_only)

        if not filename_only:
            # Use Main Product Header to get missing info
            file_path = paths[0]
            component_path = os.path.splitext(os.path.basename(file_path))[0].lower()+".xml"
            self._analyze_mph(self.read_xml_component(file_path, component_path), properties)

        return properties


class TileBasedDataProduct(BiomassBaseProduct):
    def __init__(self, product_type, zipped=None):
        pattern = [
            r"^BIO",
            product_type,
            r"(?P<mission_phase>[CIT])",
            r"G(?P<global_coverage>[_\d]{2})",
            r"T(?P<tile_number>[NS]\d{2}[EW]\d{3})",
            r"B(?P<basin_literal>[\w]{3})",
            r"(?P<baseline>\d{2})",
            r"(?P<compact_creation_date>\w{6})",
        ]
        super().__init__(product_type, filename_base_pattern=r"_".join(pattern), extension="", zipped=zipped)

    def _set_properties_from_filename(self, file_name_attrs, properties):
        # Validity start and stop cannot be recovered from filename

        biomass = properties.biomass
        biomass.type = self.product_type.split("_")[1]
        biomass.baseline = int(file_name_attrs["baseline"])
        biomass.mission_phase = file_name_attrs["mission_phase"]
        biomass.global_coverage = int(file_name_attrs["global_coverage"])
        biomass.compact_creation_date = file_name_attrs["compact_creation_date"]
        biomass.tile = file_name_attrs["tile_number"]
        biomass.basin = file_name_attrs["basin_literal"]

    def analyze(self, paths, filename_only=False):
        properties = super().analyze(paths, filename_only)

        if not filename_only:
            # Use Main Product Header to get missing info
            file_path = paths[0]
            component_path = os.path.splitext(os.path.basename(file_path))[0].lower() + ".xml"
            self._analyze_mph(self.read_xml_component(file_path, component_path), properties)

        return properties


class AuxiliaryDataFile(BiomassBaseProduct):
    def __init__(self, product_type, zipped=None):
        pattern = [
            r"^BIO",
            product_type,
            r"(?P<validity_start>\d{8}T\d{6})",
            r"(?P<validity_stop>\d{8}T\d{6})",
            r"(?P<baseline>\d{2})",
            r"(?P<compact_creation_date>\w{6})",
        ]
        super().__init__(product_type, filename_base_pattern=r"_".join(pattern), extension="", zipped=zipped)

    def _set_properties_from_filename(self, file_name_attrs, properties):
        super()._set_properties_from_filename(file_name_attrs, properties)

        biomass = properties.biomass

        biomass.type = self.product_type.split("_")[1]
        biomass.baseline = int(file_name_attrs["baseline"])
        biomass.compact_creation_date = file_name_attrs["compact_creation_date"]

    def analyze(self, paths, filename_only=False):
        properties = super().analyze(paths, filename_only)

        if not filename_only:
            # Use Main Product Header to get missing info
            file_path = paths[0]
            component_path = os.path.splitext(os.path.basename(file_path))[0].lower() + ".xml"
            self._analyze_mph(self.read_xml_component(file_path, component_path), properties)

        return properties


class BaselineConfigurationFile(BiomassBaseProduct):
    def __init__(self, product_type, zipped=None):
        pattern = [
            r"^BIO",
            r"(?P<file_class>\w{4})",
            product_type,
            r"(?P<creation_date>\d{8}T\d{6})Z",
            r"(?P<facility>\w{6})",
            r"(?P<version>\d{4})",
        ]
        super().__init__(product_type, filename_base_pattern=r"_".join(pattern), extension=".xml", zipped=zipped)

    def analyze(self, paths, filename_only=False):
        properties = super().analyze(paths, filename_only=filename_only)
        if not filename_only:
            self._analyze_eof_header(self.read_xml_component(paths[0]), properties)
        return properties

    def _set_properties_from_filename(self, file_name_attrs, properties):
        properties.core.creation_date = datetime.strptime(file_name_attrs["creation_date"], "%Y%m%dT%H%M%S")
        properties.biomass.processing_center = file_name_attrs["facility"].rstrip("_")


class Report(BiomassBaseProduct):
    def __init__(self, product_type, zipped=None):
        pattern = [
            r"^BIO",
            r"(?P<file_class>\w{4})",
            product_type,
            r"(?P<validity_start>\d{8}T\d{6})",
            r"(?P<validity_stop>\d{8}T\d{6})",
            r"(?P<originator>\w{8})",
            r"(?P<baseline>\d{2})(?P<version>\d{2})",
        ]
        super().__init__(product_type, filename_base_pattern=r"_".join(pattern), extension=".EOF", zipped=zipped)

    def _set_properties_from_filename(self, file_name_attrs, properties):
        super()._set_properties_from_filename(file_name_attrs, properties)
        properties.biomass.processing_center = file_name_attrs["originator"].rstrip("_")

    def analyze(self, paths, filename_only=False):
        properties = super().analyze(paths, filename_only)
        if not filename_only:
            self._analyze_eof_header(self.read_xml_component(paths[0]), properties)
        return properties


class QualityDisclaimerMetadataFile(BiomassBaseProduct):
    def __init__(self, product_type, zipped=None):
        pattern = [
            r"^BIO",
            r"(?P<file_class>\w{4})",
            product_type,
            r"(?P<creation_date>\d{8}T\d{6})",
            r"(?P<originator>\w{8})",
            r"(?P<identifier>\d{6})",
        ]
        super().__init__(product_type, filename_base_pattern=r"_".join(pattern), extension=".EOF", zipped=zipped)

    def analyze(self, paths, filename_only=False):
        properties = super().analyze(paths, filename_only=filename_only)
        if not filename_only:
            self._analyze_eof_header(self.read_xml_component(paths[0]), properties)
        return properties

    def _set_properties_from_filename(self, file_name_attrs, properties):
        properties.core.creation_date = datetime.strptime(file_name_attrs["creation_date"], "%Y%m%dT%H%M%S")
        properties.biomass.processing_center = file_name_attrs["originator"].rstrip("_")


class TransponderReport(BiomassBaseProduct):
    def __init__(self, product_type, zipped=None):
        self.enclose_dir = not zipped
        pattern = [
            r"^BIO",
            r"(?P<file_class>\w{4})",
            product_type,
            r"(?P<validity_start>\d{8}T\d{6})",
            r"(?P<validity_stop>\d{8}T\d{6})",
            r"(?P<version>\d{4})",
        ]
        super().__init__(product_type, filename_base_pattern=r"_".join(pattern), zipped=zipped)
        self.is_multi_file_product = True

    @property
    def use_enclosing_directory(self):
        return self.enclose_dir

    def enclosing_directory(self, properties):
        return properties.core.product_name

    def identify(self, paths):
        for path in paths:
            if os.path.isdir(path):
                return False
            if re.match(self.filename_pattern, os.path.basename(path)) is None:
                return False
        return True

    def analyze(self, paths, filename_only=False):
        if len(paths) > 1:
            # make sure we pass the header file as first entry for analyze()
            header_files = [path for path in paths if path.endswith(".HDR")]
            if len(header_files) == 0:
                raise muninnError("product does not contain a .HDR file")
            header_file = header_files[0]
            paths.remove(header_file)
            paths.insert(0, header_file)

        properties = super().analyze(paths, filename_only)

        if not filename_only:
            # Use header file to extract info
            if len(paths) == 1 and self.is_zipped(paths[0]):
                component_path = os.path.splitext(os.path.basename(paths[0]))[0] + ".HDR"
            else:
                component_path = None
            self._analyze_eof_header(self.read_xml_component(paths[0], component_path), properties)

        return properties

    def _set_properties_from_filename(self, file_name_attrs, properties):
        super()._set_properties_from_filename(file_name_attrs, properties)

        properties.biomass.type = self.product_type[4:7]


class EOFFile(BiomassBaseProduct):
    def __init__(self, product_type, filename_base_pattern=None, zipped=None):
        if filename_base_pattern is None:
            pattern = [
                r"^BIO",
                r"(?P<file_class>\w{4})",
                product_type,
                r"(?P<validity_start>\d{8}T\d{6})",
                r"(?P<validity_stop>\d{8}T\d{6})",
                r"(?P<version>\d{4})",
            ]
            filename_base_pattern = r"_".join(pattern)
        super().__init__(product_type, filename_base_pattern=filename_base_pattern, extension=".EOF", zipped=zipped)

    def analyze(self, paths, filename_only=False):
        properties = super().analyze(paths, filename_only)

        if not filename_only:
            # Use content of Fixed_Header
            self._analyze_eof_header(self.read_xml_component(paths[0]), properties)

        return properties


class TransponderEOFFile(EOFFile):
    """Only difference with EOFFile is that biomass.type is set """
    def _set_properties_from_filename(self, file_name_attrs, properties):
        super()._set_properties_from_filename(file_name_attrs, properties)
        properties.biomass.type = self.product_type[4:7]


class EOFFileWithBaseline(EOFFile):
    """Only difference with EOFFile is that baseline is a component of the version (which is always 0)"""
    def __init__(self, product_type, zipped=None):
        pattern = [
            r"^BIO",
            r"(?P<file_class>\w{4})",
            product_type,
            r"(?P<validity_start>\d{8}T\d{6})",
            r"(?P<validity_stop>\d{8}T\d{6})",
            r"(?P<baseline>00)(?P<version>\d{2})",
        ]
        super().__init__(product_type, filename_base_pattern=r"_".join(pattern), zipped=zipped)


class EOFFileWithAcquisition(EOFFile):
    """Difference with EOFFileWithBaseline is the addition of acquisition_station"""
    def __init__(self, product_type, zipped=None):
        pattern = [
            r"^BIO",
            r"(?P<file_class>\w{4})",
            product_type,
            r"(?P<validity_start>\d{8}T\d{6})",
            r"(?P<validity_stop>\d{8}T\d{6})",
            r"(?P<acquisition_station>\w{8})",
            r"(?P<baseline>00)(?P<version>\d{2})",
        ]
        super().__init__(product_type, filename_base_pattern=r"_".join(pattern), zipped=zipped)

    def _set_properties_from_filename(self, file_name_attrs, properties):
        super()._set_properties_from_filename(file_name_attrs, properties)
        properties.biomass.acquisition_station = file_name_attrs["acquisition_station"]


class EOFFileWithPhase(EOFFile):
    """Difference with EOFFile is the addition of mission_phase"""
    def __init__(self, product_type, zipped=None):
        pattern = [
            r"^BIO",
            r"(?P<file_class>\w{4})",
            product_type,
            r"(?P<validity_start>\d{8}T\d{6})",
            r"(?P<validity_stop>\d{8}T\d{6})",
            r"(?P<mission_phase>\w{8})",
            r"(?P<version>\d{4})",
        ]
        super().__init__(product_type, filename_base_pattern=r"_".join(pattern), zipped=zipped)

    def _set_properties_from_filename(self, file_name_attrs, properties):
        super()._set_properties_from_filename(file_name_attrs, properties)
        properties.biomass.mission_phase = file_name_attrs["mission_phase"][0]


_product_types = dict(
    [(product_type, FrameBasedDataProduct(product_type)) for product_type in L0_PRODUCT_TYPES] +
    [(product_type, FrameBasedDataProduct(product_type)) for product_type in L1_PRODUCT_TYPES] +
    [(product_type, FrameBasedDataProduct(product_type)) for product_type in L2A_PRODUCT_TYPES] +
    [(product_type, TileBasedDataProduct(product_type)) for product_type in L2B_L3_PRODUCT_TYPES] +
    [(product_type, AuxiliaryDataFile(product_type)) for product_type in ADF_PRODUCT_TYPES] +
    [(product_type, BaselineConfigurationFile(product_type)) for product_type in BASELINE_CFG_PRODUCT_TYPES] +
    [(product_type, Report(product_type)) for product_type in REPORT_PRODUCT_TYPES] +
    [(product_type, QualityDisclaimerMetadataFile(product_type)) for product_type in QUALITY_DISCLAIMER_PRODUCT_TYPES] +
    [(product_type, TransponderEOFFile(product_type)) for product_type in TRANSPONDER_EOF_FILE_TYPES] +
    [(product_type, TransponderReport(product_type)) for product_type in TRANSPONDER_REPORT_TYPES] +
    [(product_type, EOFFile(product_type)) for product_type in MISSION_PLANNING_EOF_TYPES] +
    [(product_type, EOFFileWithBaseline(product_type)) for product_type in MISSION_PLANNING_EOF_WITH_BASELINE_TYPES] +
    [(product_type, EOFFileWithAcquisition(product_type)) for product_type in MISSION_PLANNING_EOF_WITH_ACQ_TYPES] +
    [(product_type, EOFFileWithPhase(product_type)) for product_type in MISSION_PLANNING_EOF_WITH_PHASE_TYPES]
)


def product_types():
    return _product_types.keys()


def product_type_plugin(product_type):
    return _product_types.get(product_type)
