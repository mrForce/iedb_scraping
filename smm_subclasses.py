#!/usr/bin/env python

#
# Generated Wed May  3 10:52:22 2017 by generateDS.py version 2.25a.
#
# Command line options:
#   ('-o', 'smm.py')
#   ('-s', 'smm_subclasses.py')
#   ('--super', 'smm_api')
#
# Command line arguments:
#   smm.xsd
#
# Command line:
#   /usr/bin/generateDS.py -o "smm.py" -s "smm_subclasses.py" --super="smm_api" smm.xsd
#
# Current working directory (os.getcwd()):
#   smmlinux
#

import sys
from lxml import etree as etree_

import smm_api as supermod

def parsexml_(infile, parser=None, **kwargs):
    if parser is None:
        # Use the lxml ElementTree compatible parser so that, e.g.,
        #   we ignore comments.
        parser = etree_.ETCompatXMLParser()
    doc = etree_.parse(infile, parser=parser, **kwargs)
    return doc

#
# Globals
#

ExternalEncoding = 'ascii'

#
# Data representation classes
#


class SMMTrainInputSub(supermod.SMMTrainInput):
    def __init__(self, OutputFile=None, TrainingData=None, MatrixCalculation=None, PairCalculation=None):
        super(SMMTrainInputSub, self).__init__(OutputFile, TrainingData, MatrixCalculation, PairCalculation, )
supermod.SMMTrainInput.subclass = SMMTrainInputSub
# end class SMMTrainInputSub


class SMMTrainOutputSub(supermod.SMMTrainOutput):
    def __init__(self, Log=None, SMMPredictor=None):
        super(SMMTrainOutputSub, self).__init__(Log, SMMPredictor, )
supermod.SMMTrainOutput.subclass = SMMTrainOutputSub
# end class SMMTrainOutputSub


class SMMParametersSub(supermod.SMMParameters):
    def __init__(self, Repeats=None, LambdaRange=None, Precission=None, SeedRandomizer=None, MaxNormIterations=None):
        super(SMMParametersSub, self).__init__(Repeats, LambdaRange, Precission, SeedRandomizer, MaxNormIterations, )
supermod.SMMParameters.subclass = SMMParametersSub
# end class SMMParametersSub


class SMMPredictorSub(supermod.SMMPredictor):
    def __init__(self, Alphabet=None, SequenceLength=None, SeqMatrix=None, SeqPair=None):
        super(SMMPredictorSub, self).__init__(Alphabet, SequenceLength, SeqMatrix, SeqPair, )
supermod.SMMPredictor.subclass = SMMPredictorSub
# end class SMMPredictorSub


class SMMPredictInputSub(supermod.SMMPredictInput):
    def __init__(self, OutputFile=None, SMMPredictor=None, Predict=None):
        super(SMMPredictInputSub, self).__init__(OutputFile, SMMPredictor, Predict, )
supermod.SMMPredictInput.subclass = SMMPredictInputSub
# end class SMMPredictInputSub


class SMMPredictOutputSub(supermod.SMMPredictOutput):
    def __init__(self, Log=None, Predict=None):
        super(SMMPredictOutputSub, self).__init__(Log, Predict, )
supermod.SMMPredictOutput.subclass = SMMPredictOutputSub
# end class SMMPredictOutputSub


class SeqMatrixSub(supermod.SeqMatrix):
    def __init__(self, Offset=None, MatCoef=None):
        super(SeqMatrixSub, self).__init__(Offset, MatCoef, )
supermod.SeqMatrix.subclass = SeqMatrixSub
# end class SeqMatrixSub


class LogSub(supermod.Log):
    def __init__(self, success=None, LogText=None, ErrorMessage=None):
        super(LogSub, self).__init__(success, LogText, ErrorMessage, )
supermod.Log.subclass = LogSub
# end class LogSub


class TrainingDataTypeSub(supermod.TrainingDataType):
    def __init__(self, Alphabet=None, SequenceLength=None, SequenceData=None, SeqMatrix=None):
        super(TrainingDataTypeSub, self).__init__(Alphabet, SequenceLength, SequenceData, SeqMatrix, )
supermod.TrainingDataType.subclass = TrainingDataTypeSub
# end class TrainingDataTypeSub


class SequenceDataTypeSub(supermod.SequenceDataType):
    def __init__(self, DataPoint=None):
        super(SequenceDataTypeSub, self).__init__(DataPoint, )
supermod.SequenceDataType.subclass = SequenceDataTypeSub
# end class SequenceDataTypeSub


class DataPointTypeSub(supermod.DataPointType):
    def __init__(self, Sequence=None, Threshold=None, Measured=None):
        super(DataPointTypeSub, self).__init__(Sequence, Threshold, Measured, )
supermod.DataPointType.subclass = DataPointTypeSub
# end class DataPointTypeSub


class ThresholdTypeSub(supermod.ThresholdType):
    def __init__(self, Greater=None, Lesser=None):
        super(ThresholdTypeSub, self).__init__(Greater, Lesser, )
supermod.ThresholdType.subclass = ThresholdTypeSub
# end class ThresholdTypeSub


class MatrixCalculationTypeSub(supermod.MatrixCalculationType):
    def __init__(self, UseExternalMatrix=None, SMMMatrix=None):
        super(MatrixCalculationTypeSub, self).__init__(UseExternalMatrix, SMMMatrix, )
supermod.MatrixCalculationType.subclass = MatrixCalculationTypeSub
# end class MatrixCalculationTypeSub


class UseExternalMatrixTypeSub(supermod.UseExternalMatrixType):
    def __init__(self, AdjustToSequenceData=None, AdjustOffsetToSequenceData=None):
        super(UseExternalMatrixTypeSub, self).__init__(AdjustToSequenceData, AdjustOffsetToSequenceData, )
supermod.UseExternalMatrixType.subclass = UseExternalMatrixTypeSub
# end class UseExternalMatrixTypeSub


class SMMMatrixTypeSub(supermod.SMMMatrixType):
    def __init__(self, LambdaGrouping=None, SMMParameters=None):
        super(SMMMatrixTypeSub, self).__init__(LambdaGrouping, SMMParameters, )
supermod.SMMMatrixType.subclass = SMMMatrixTypeSub
# end class SMMMatrixTypeSub


class LambdaGroupingTypeSub(supermod.LambdaGroupingType):
    def __init__(self, Individual=None, Position=None, One=None):
        super(LambdaGroupingTypeSub, self).__init__(Individual, Position, One, )
supermod.LambdaGroupingType.subclass = LambdaGroupingTypeSub
# end class LambdaGroupingTypeSub


class PairCalculationTypeSub(supermod.PairCalculationType):
    def __init__(self, IndividualLambdas=None, PairCoefCriteria=None, SMMParameters=None):
        super(PairCalculationTypeSub, self).__init__(IndividualLambdas, PairCoefCriteria, SMMParameters, )
supermod.PairCalculationType.subclass = PairCalculationTypeSub
# end class PairCalculationTypeSub


class PairCoefCriteriaTypeSub(supermod.PairCoefCriteriaType):
    def __init__(self, MinCount=None, MaxDisagrement=None):
        super(PairCoefCriteriaTypeSub, self).__init__(MinCount, MaxDisagrement, )
supermod.PairCoefCriteriaType.subclass = PairCoefCriteriaTypeSub
# end class PairCoefCriteriaTypeSub


class RepeatsTypeSub(supermod.RepeatsType):
    def __init__(self, Bagging=None, CrossValidation=None):
        super(RepeatsTypeSub, self).__init__(Bagging, CrossValidation, )
supermod.RepeatsType.subclass = RepeatsTypeSub
# end class RepeatsTypeSub


class LambdaRangeTypeSub(supermod.LambdaRangeType):
    def __init__(self, Min=None, Start=None, Max=None):
        super(LambdaRangeTypeSub, self).__init__(Min, Start, Max, )
supermod.LambdaRangeType.subclass = LambdaRangeTypeSub
# end class LambdaRangeTypeSub


class SeqPairTypeSub(supermod.SeqPairType):
    def __init__(self, PairCoef=None):
        super(SeqPairTypeSub, self).__init__(PairCoef, )
supermod.SeqPairType.subclass = SeqPairTypeSub
# end class SeqPairTypeSub


class PairCoefTypeSub(supermod.PairCoefType):
    def __init__(self, Letter1=None, Position1=None, Letter2=None, Position2=None, Value=None):
        super(PairCoefTypeSub, self).__init__(Letter1, Position1, Letter2, Position2, Value, )
supermod.PairCoefType.subclass = PairCoefTypeSub
# end class PairCoefTypeSub


class PredictTypeSub(supermod.PredictType):
    def __init__(self, Sequence=None):
        super(PredictTypeSub, self).__init__(Sequence, )
supermod.PredictType.subclass = PredictTypeSub
# end class PredictTypeSub


class PredictType1Sub(supermod.PredictType1):
    def __init__(self, ScanLength=None, Sequence=None, Predictions=None):
        super(PredictType1Sub, self).__init__(ScanLength, Sequence, Predictions, )
supermod.PredictType1.subclass = PredictType1Sub
# end class PredictType1Sub


class MatCoefTypeSub(supermod.MatCoefType):
    def __init__(self, Letter=None, Position=None, Value=None):
        super(MatCoefTypeSub, self).__init__(Letter, Position, Value, )
supermod.MatCoefType.subclass = MatCoefTypeSub
# end class MatCoefTypeSub


def get_root_tag(node):
    tag = supermod.Tag_pattern_.match(node.tag).groups()[-1]
    rootClass = None
    rootClass = supermod.GDSClassesMapping.get(tag)
    if rootClass is None and hasattr(supermod, tag):
        rootClass = getattr(supermod, tag)
    return tag, rootClass


def parse(inFilename, silence=False):
    parser = None
    doc = parsexml_(inFilename, parser)
    rootNode = doc.getroot()
    rootTag, rootClass = get_root_tag(rootNode)
    if rootClass is None:
        rootTag = 'SMMTrainInput'
        rootClass = supermod.SMMTrainInput
    rootObj = rootClass.factory()
    rootObj.build(rootNode)
    # Enable Python to collect the space used by the DOM.
    doc = None
    if not silence:
        sys.stdout.write('<?xml version="1.0" ?>\n')
        rootObj.export(
            sys.stdout, 0, name_=rootTag,
            namespacedef_='',
            pretty_print=True)
    return rootObj


def parseEtree(inFilename, silence=False):
    parser = None
    doc = parsexml_(inFilename, parser)
    rootNode = doc.getroot()
    rootTag, rootClass = get_root_tag(rootNode)
    if rootClass is None:
        rootTag = 'SMMTrainInput'
        rootClass = supermod.SMMTrainInput
    rootObj = rootClass.factory()
    rootObj.build(rootNode)
    # Enable Python to collect the space used by the DOM.
    doc = None
    mapping = {}
    rootElement = rootObj.to_etree(None, name_=rootTag, mapping_=mapping)
    reverse_mapping = rootObj.gds_reverse_node_mapping(mapping)
    if not silence:
        content = etree_.tostring(
            rootElement, pretty_print=True,
            xml_declaration=True, encoding="utf-8")
        sys.stdout.write(content)
        sys.stdout.write('\n')
    return rootObj, rootElement, mapping, reverse_mapping


def parseString(inString, silence=False):
    from StringIO import StringIO
    parser = None
    doc = parsexml_(StringIO(inString), parser)
    rootNode = doc.getroot()
    rootTag, rootClass = get_root_tag(rootNode)
    if rootClass is None:
        rootTag = 'SMMTrainInput'
        rootClass = supermod.SMMTrainInput
    rootObj = rootClass.factory()
    rootObj.build(rootNode)
    # Enable Python to collect the space used by the DOM.
    doc = None
    if not silence:
        sys.stdout.write('<?xml version="1.0" ?>\n')
        rootObj.export(
            sys.stdout, 0, name_=rootTag,
            namespacedef_='')
    return rootObj


def parseLiteral(inFilename, silence=False):
    parser = None
    doc = parsexml_(inFilename, parser)
    rootNode = doc.getroot()
    rootTag, rootClass = get_root_tag(rootNode)
    if rootClass is None:
        rootTag = 'SMMTrainInput'
        rootClass = supermod.SMMTrainInput
    rootObj = rootClass.factory()
    rootObj.build(rootNode)
    # Enable Python to collect the space used by the DOM.
    doc = None
    if not silence:
        sys.stdout.write('#from smm_api import *\n\n')
        sys.stdout.write('import smm_api as model_\n\n')
        sys.stdout.write('rootObj = model_.rootClass(\n')
        rootObj.exportLiteral(sys.stdout, 0, name_=rootTag)
        sys.stdout.write(')\n')
    return rootObj


USAGE_TEXT = """
Usage: python ???.py <infilename>
"""


def usage():
    print(USAGE_TEXT)
    sys.exit(1)


def main():
    args = sys.argv[1:]
    if len(args) != 1:
        usage()
    infilename = args[0]
    parse(infilename)


if __name__ == '__main__':
    #import pdb; pdb.set_trace()
    main()
