import numpy as np
import copy
from scipy.io import mmread, mmwrite
from scipy.sparse import coo_matrix

def make_realm_from_complexm(infile, outfile, OUT_STD=False):
    # Read file of Matrix Market format, and get rank of matrix
    # (Asume rank of matrix is in the second line)
    if OUT_STD:
        print("Reading file {}...".format(infile))
        print()

    with open(infile) as f:
        lines = f.readlines()

    nums = [int(x) for x in lines[1].split()]
    ranka = nums[0]

    if OUT_STD:
        print("Rank of input matrix is {}.".format(ranka))
        print()

    # Read file of Matrix Market format
    a = mmread(infile)
    if OUT_STD:
        print("The input matrix is:")
        print(type(a))
        print(a)
        print()

    # Make output real matrix from complex matrix
    if OUT_STD:
        print("Rank of output matrix is {}.".format(ranka*2))
        print()

    a_real_data = [np.real(x) for x in a.data]
    a_imag_data = [np.imag(x) for x in a.data]
    #print(a_real_data)
    #print(a_imag_data)

    # P + i Q -> | P(1) -Q(3) |
    #            | Q(2)  P(4) |

    #(1)
    outdata_tmp = copy.copy(a_real_data)
    outrow_tmp = list(a.row)
    outcol_tmp = list(a.col)

    #(2)
    tmpdata = copy.copy(a_imag_data)
    tmprow = [x+ranka for x in a.row]
    tmpcol = list(a.col)
    outdata_tmp.extend(tmpdata)
    outrow_tmp.extend(tmprow)
    outcol_tmp.extend(tmpcol)

    #(3)
    tmpdata = [-x for x in a_imag_data]
    tmprow = list(a.row)
    tmpcol = [x+ranka for x in a.col]
    outdata_tmp.extend(tmpdata)
    outrow_tmp.extend(tmprow)
    outcol_tmp.extend(tmpcol)
    
    #(4)
    tmpdata = copy.copy(a_real_data)
    tmprow = [x+ranka for x in a.row]
    tmpcol = [x+ranka for x in a.col]
    outdata_tmp.extend(tmpdata)
    outrow_tmp.extend(tmprow)
    outcol_tmp.extend(tmpcol)

    outdata = []
    outrow = []
    outcol = []
    for i, x in enumerate(outdata_tmp):
        if x != 0.0:
            outdata.append(x)
            outrow.append(outrow_tmp[i])
            outcol.append(outcol_tmp[i])
   
    b = coo_matrix((outdata, (outrow, outcol)), shape=(ranka*2, ranka*2))
    if OUT_STD:
        print("The output matrix is:")
        print(type(b))
        print(b)
        print()
        print("Output to file {}...".format(outfile))
    mmwrite(outfile, b)

if __name__ == '__main__':
    infile = "test.mtx"
    outfile = "output.mtx"
    make_realm_from_complexm(infile, outfile, OUT_STD=True)
