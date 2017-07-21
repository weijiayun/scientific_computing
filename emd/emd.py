from IMF import SoloIMF, BiIMF
import numpy as np
from scipy.signal import hilbert

def hilbertSpectrum(signal, fs):
    hilbertSignal = hilbert(signal)
    instantaneous_phase = np.unwrap(np.angle(hilbertSignal))

    instantaneous_frequency = (np.diff(instantaneous_phase) /
                                   (2.0 * np.pi) * fs)

    return instantaneous_frequency, hilbertSignal

def emd(signal, numOfIMF, fittingType='cs',  SDLowerLimt=0, SDUpperLimit=0.3):
    soloIMFInst = SoloIMF(signal, numOfIMF, fittingType=fittingType, SDLowerLimt=SDLowerLimt, SDUpperLimit=SDUpperLimit)
    return soloIMFInst.getIMF()

def bi_emd(signal, mComponentH=1, nComponentV=1, fittingType='cs',  SDLowerLimt=0, SDUpperLimit=0.3):
    biIMFInst = BiIMF(signal, mComponentH, nComponentV, fittingType=fittingType, SDUpperLimit=SDUpperLimit, SDLowerLimt=SDLowerLimt)
    return biIMFInst.getIMF()


def eemd(signal, numOfIMF, numOfEnsemble=None, ensembleRadio=0.2, fittingType='cs',  SDLowerLimt=0, SDUpperLimit=0.3):

    signalSize = np.max(signal.shape)
    if not numOfEnsemble:
        numOfEnsemble = int(np.ceil(np.log2(signalSize)))

    ystd = np.std(signal)
    std_signal = signal / ystd

    tmpEnsemble = np.zeros((numOfEnsemble, numOfIMF, signalSize), dtype='float64')

    for iii in range(numOfEnsemble):
        randomArr = np.random.standard_normal((signalSize,))

        signal1 = std_signal + randomArr * ensembleRadio
        signal2 = std_signal - randomArr * ensembleRadio

        soloIMFInst1 = SoloIMF(signal1, numOfIMF, fittingType=fittingType, SDLowerLimt=SDLowerLimt, SDUpperLimit=SDUpperLimit)
        soloIMFInst2 = SoloIMF(signal2, numOfIMF, fittingType=fittingType, SDLowerLimt=SDLowerLimt, SDUpperLimit=SDUpperLimit)

        imfs1 = soloIMFInst1.getIMF()
        imfs2 = soloIMFInst2.getIMF()

        for i in range(numOfIMF):
            tmpEnsemble[iii, i, :] = (imfs1[i].getCi() + imfs2[i].getCi())* ystd/2

    container = np.zeros((numOfIMF, signalSize))
    for i in range(numOfIMF):
        container[i, :] = np.mean(tmpEnsemble[:, i, :], axis=0)

    return container

def bi_eemd(signal, mComponentH=1, nComponentV=1, numOfEnsemble=20, ensembleRadio=0.2, fittingType='cs',  SDLowerLimt=0, SDUpperLimit=0.3):

    signalSize = np.max(signal.shape)
    if not numOfEnsemble:
        numOfEnsemble = int(np.ceil(np.log2(signalSize)))

    numOfIMF = np.min([mComponentH, nComponentV])

    ystd = np.std(signal)

    std_signal = signal / ystd

    tmpEnsemble = np.zeros((numOfEnsemble, numOfIMF, signal.shape[0], signal.shape[1]), dtype='float64')

    for iii in range(numOfEnsemble):
        randomArr = np.random.standard_normal(signal.shape)
        signal1 = std_signal + randomArr * ensembleRadio
        signal2 = std_signal - randomArr * ensembleRadio

        soloIMFInst1 = BiIMF(signal1, mComponentH, nComponentV, fittingType=fittingType, SDLowerLimt=SDLowerLimt, SDUpperLimit=SDUpperLimit)
        soloIMFInst2 = BiIMF(signal2, mComponentH, nComponentV, fittingType=fittingType, SDLowerLimt=SDLowerLimt, SDUpperLimit=SDUpperLimit)

        imfs1 = soloIMFInst1.getIMF()
        imfs2 = soloIMFInst2.getIMF()

        for i in range(numOfIMF):
            tmpEnsemble[iii, i, :, :] = (imfs1[i, :, :] + imfs2[i, :, :])* ystd/2

    container = np.zeros((numOfIMF, signal.shape[0], signal.shape[1]))
    for i in range(numOfIMF):
        container[i, :, :] = np.mean(tmpEnsemble[:, i, :, :], axis=0)

    return container






