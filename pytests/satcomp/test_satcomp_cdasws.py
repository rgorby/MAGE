#test_satcomp_cdasws
import datetime
import kaipy.satcomp.scutils as scutils
import pytest
import numpy as np

def test_getscIds():
    scIdDict = scutils.getScIds()
    assert type(scIdDict) == dict, "Returned type is {}, but should be type dict".format(type(scIdDict))

@pytest.mark.parametrize("key",[key for key in scutils.getScIds().keys()])
def test_emphem(key):

    scIdDict = scutils.getScIds()
    assert 'Ephem' in scIdDict[key], \
        '{} lacks Ephem key'.format(key)
    assert 'Id' in scIdDict[key]['Ephem'], \
        '{} lacks Ephem Id key'.format(key)
    assert 'Data' in scIdDict[key]['Ephem'], \
        '{} lacks Ephem Data key'.format(key)
    assert 'CoordSys' in scIdDict[key]['Ephem'], \
        '{} lacks Ephem CoordSys key'.format(key)


@pytest.mark.parametrize("key",[key for key in scutils.getScIds().keys()])
def test_pullVar(key):

    scIdDict = scutils.getScIds()

    for var in scIdDict[key]:
        if var != "_testing":
            tStart, tEnd = scutils.getCdasDsetInterval(scIdDict[key][var]['Id'])
            assert tStart is not None,\
                '{} did not have valid start time'.format(key)
            # t0 = tStart
            # t0dt = datetime.datetime.strptime(t0, "%Y-%m-%dT%H:%M:%S.%fZ")
            # t1 = (t0dt + datetime.timedelta(days=1)).strftime("%Y-%m-%dT%H:%M:%S.%fZ")
            t1 = tEnd
            t1dt = datetime.datetime.strptime(t1, "%Y-%m-%dT%H:%M:%S.%fZ")
            t0 = (t1dt - datetime.timedelta(days=1)).strftime("%Y-%m-%dT%H:%M:%S.%fZ")
            dset_id = scIdDict[key][var]['Id']
            dset_vname = scIdDict[key][var]['Data']
            status,data = scutils.pullVar(scIdDict[key][var]['Id'],scIdDict[key][var]['Data'],
                      t0,t1,60.0)
            assert status['http']['status_code'] == 200, \
                "pullVar failed to return for {},{}".format(key,var)
            if isinstance(scIdDict[key][var]['Data'],str):
                assert scIdDict[key][var]['Data'] in data.keys(), \
                    "pullVar error - data does not contain {} for {},{}".format(scIdDict[key][var]['Data'],key,var)
            if isinstance(scIdDict[key][var]['Data'],list):
                for item in scIdDict[key][var]['Data']:
                    assert item in data.keys(), \
                        "pullVar error - data does not contain {} for {},{}".format(item,key,var)

                    
            

@pytest.fixture
def exampleObs():
    return np.array([3, -0.5, 2, 7])
@pytest.fixture
def examplePred():
    return np.array([2.5, 0.0, 2, 8])

def test_computeErrors(exampleObs,examplePred):
    MAE,MSE,RMSE,MAPE,RSE,PE = scutils.computeErrors(exampleObs,examplePred)
    assert MAE == pytest.approx(0.5)
    assert MSE == pytest.approx(0.375)
    assert RMSE == pytest.approx(0.6123724356957945)
    assert MAPE == pytest.approx(0.3273809523809524)
    assert RSE == pytest.approx(0.05139186295503212)
    assert PE == pytest.approx(0.9486081370449679)

    return
