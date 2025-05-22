
from PhotDG import PhotDG
from SpecDG import SpecDG

if __name__ == '__main__':
    # get the DG
    dg = PhotDG()

    # print the DG parameters
    dg.print_parameters()

    # get the photometry for a grid point
    dg.photGet('RV31_BC6-WD01', 'dusty', 'burst', 0.004, 10.0, 1e6, 1.0)

    # now for the spectrum
    sg = SpecDG()
    gid = sg.findGidFromParam('RV31_BC6-WD01', 'dusty', 'burst', 0.004, 10.0, 1e6, 1.0)
    filename = sg.findFile(gid)
    print(gid, filename)
    sg.specGet(filename)
    sg.specPlot()

    # plot that SED
    dg.photPlot()
