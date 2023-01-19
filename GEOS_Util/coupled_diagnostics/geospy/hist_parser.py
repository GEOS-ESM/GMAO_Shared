import re

def get_history(fname, GRID_LABELS=True):
    '''
    Parses HIHSTORY.rc file given by 'fname' and returns a
    python dictionary which contains a header with history options
    and a list of collections.
    '''
    
    history=dict(header={}, collections=[])
    
    with open(fname) as ff:
        # Read file, delete comments
        intext=re.sub('#.*','',ff.read())
    
    # Split input into '::' blocks, remove leading and trailing blanks
    blocks=[s.strip() for s in re.split('\W::\W',intext)]
    # Remove empty blocks 
    blocks=list(filter(lambda x: x!='', blocks))
    
    # Process header
    headstr=blocks[0]
    header={k:v for (k,v) in re.findall('(.*):(.*)', headstr)}
    
    # Process collection list
    colstr=re.search('COLLECTIONS:(.*)',headstr, re.DOTALL).group(1)
    collist=[s.strip(" \t'") for s in colstr.split('\n')]
    collist=list(filter(lambda x: x!='', collist)) # Remove blank entries
    header.update(COLLECTIONS=collist)
    history.update(header=header)
    
    # Process collections
    start_block = 3 if GRID_LABELS else 1
    for colstr in blocks[start_block:]:
        colname=re.match('.*?\.',colstr).group().strip('.')
        print(colname)
        if colname in header['COLLECTIONS']:
            colopts={k.replace(f'{colname}.','').strip():v.strip(" \t,'") for (k,v) in re.findall('(.*):(.*)',colstr)}
            
            # Process var list
            varstr=re.search('fields:(.*)',colstr, re.DOTALL).group(1)
            varlist=[s.strip(" \t,") for s in varstr.split('\n')]
            varlist=list(filter(lambda x: x!='', varlist)) # remove blank entries
            colopts.update(name=colname, fields=varlist)
            history['collections'].append(colopts)
        else:
            print(f'{colname} is not in collection list')
            
    return history
