import futils

def g2g(datain,tiledata,hdimsout):
    '''
    Interpolates field from one horizontal grid onto another using exchange 
    grid (tiles).
    
    datain - input array (at least 2d)
    tiledata - array of records of size Nt (Nt - number of tiles).
    tiledata['iin'], tiledata['jin'], tiledata['iout'], 
    tiledata['jout'], tiledata['frac']  - 
    arrays of size Nt of inexes on the input grid, indexes on the output grid, 
    corresponding to a tile N, and fraction of tile on the output grid
    hdimsout - (imout,jmout), dimensions of output grid
    '''
    shin=datain.shape
    shout=datain.shape[:-2]+hdimsout
            
    return sp.ma.masked_values(
        futils.g2g(datain.filled().reshape((-1,)+shin[-2:]),
                   tiledata['iin'].astype(int),
                   tiledata['jin'].astype(int),
                   tiledata['iout'].astype(int),
                   tiledata['jout'].astype(int),
                   tiledata['frac'],
                   datain.fill_value,
                   shout[-1],
                   shout[-2]),datain.fill_value).reshape(shout)

def g2g_field(field,ogrid,tiledata):
    '''
    Interpolates field from one horizontal grid onto another using 
    exchange grid (tiles).
    
    ogrid - output grid

    tiledata - array of records of size Nt (Nt - number of tiles).
    tiledata['iin'], tiledata['jin'], tiledata['iout'], tiledata['jout'], 
    tiledata['frac']  - 
    arrays of size Nt of indexes on the input grid, indexes on the output grid, 
    corresponding to a tile N, and fraction of tile on the output grid.
    '''
    outdims=ogrid.dims[-2:]

    field.data=utl.g2g(field.data,tiledata,outdims)
    field.grid=grid.Grid(ogrid['lon'],ogrid['lat'],field.grid['lev'])

