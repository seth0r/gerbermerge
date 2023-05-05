#!/usr/bin/env python3
import os
import sys
import re
import glob

class Gerber_RS_274X():
    def __init__(self, fname):
        self.fname = fname
        self.extra_params = {}
        self.ap_macros = {}
        self.ap_descs = {}
        self.layers = []
        self.bbox = None

        self.units = None
        self.z_tail_trim = False
        self.f_x = None
        self.f_y = None

        self.read()
        self.calc_bbox()

        if self.units and self.f_x and self.f_y:
            if len(self.ap_descs) + len(self.layers) > 0:
                print("%s: OK %s" % (fname,self.bbox))
            else:
                print("%s: empty" % fname)
        else:
            print(self)
            raise Exception("Parameters missing")

    def __str__(self):
        return "{1}:  Units: {0.units}, Tailing zeroes trimmed: {0.z_tail_trim}, Format: x{0.f_x} y{0.f_y}, AP macros: {2}, AP descs: {3}".format(
            self, self.fname, len(self.ap_macros), len(self.ap_descs) ) + "\n" + str(self.extra_params)

    def write(self,f):
        self.write_head(f)
        f.write("\n\n")
        self.write_layers(f)
        f.write("M02*")

    def write_head(self,f):
        f.write("G04 GerberMerger*\n")
        f.write("G75*\n")
        f.write("%%MO%s*%%\n" % self.units.upper())
        f.write("%%FSLAX%d%dY%d%d*%%\n" % ( self.f_x[0], self.f_x[1], self.f_y[0], self.f_y[1] ))
        for c,pl in self.extra_params.items():
            f.write("%%%s%s*%%\n" % ( c,"*".join(pl) ))
        for n,m in self.ap_macros.items():
            f.write("%%AM%s*%s*%%\n" % ( n,m ))
        f.write("G01*\n")
        for n,d in self.ap_descs.items():
            f.write("%%ADD%s%s*%%\n" % ( n,d ))

    def write_layers(self,f):
        for c,l in self.layers:
            f.write("%s*\n" % c)
            for x,y,d in l:
                sx = str(int(x * (10**self.f_x[1])))
                sy = str(int(y * (10**self.f_y[1])))
                f.write("X%sY%s%s*\n" % ( sx,sy,d ))

    def read(self):
        f = open(self.fname,"rt")
        l = f.readline()
        while l:
            l = l.strip()
            if l.startswith("%"):
                while l.count("%") % 2 > 0:
                    l += f.readline().strip()
                for p in filter(lambda e: e, l.split("%")):
                    cmd = p[:2].upper()
                    blocks = filter(lambda e: e, p[2:].split("*"))
                    if hasattr(self,"x_%s" % cmd):
                        getattr(self,"x_%s" % cmd)( *blocks )
                    else:
                        self.extra_params[ cmd ] = list(blocks)
            else:
                for code in filter(lambda e: e, l.upper().split("*")):
                    c = code[:3]
                    if c[0] == 'G':
                        if len(self.ap_descs) == 0 and c in ["G04","G75","G01"]:
                            continue
                        if hasattr(self,"g_%s" % code):
                            getattr(self,"g_%s" % code)()
                        else:
                            raise Exception("G code %s unknown." % c)
                    elif c[0] == 'M':
                        if hasattr(self,"m_%s" % code):
                            getattr(self,"m_%s" % code)()
                        else:
                            raise Exception("M code %s unknown." % code)
                    elif c[0] == 'D':
                        if int(code[1:]) in self.ap_descs:
                            self.layers.append( (code,[]) )
                        else:
                            raise Exception("Aperture %s undefined" % code)
                    elif c[0] == 'N':
                        pass
                    elif c[0] == 'X':
                        self.add_coord( code )
                    else:
                        raise Exception("%s unknown." % code)
            l = f.readline()
        f.close()

    def add_coord(self, coord):
        for c in coord:
            if c not in "XYD0123456789-":
                raise Exception("Unknown coordinate:", coord)
        coord,_,d = coord.rpartition("D")
        x,_,y = coord.partition("Y")
        x = x.strip("X")
        if self.z_tail_trim:
            Exception("Format without tailing zeroes not implemented")
        else:
            x = int(x) / (10**self.f_x[1])
            y = int(y) / (10**self.f_y[1])
        self.layers[-1][1].append( ( x,y,"D%s" % d ) )

    def calc_bbox(self):
        bbox = [None,None,None,None]
        for c,l in self.layers:
            for x,y,d in l:
                if bbox[0] is None or x < bbox[0]:
                    bbox[0] = x
                if bbox[2] is None or x > bbox[2]:
                    bbox[2] = x
                if bbox[1] is None or y < bbox[1]:
                    bbox[1] = y
                if bbox[3] is None or y > bbox[3]:
                    bbox[3] = y
        if None not in bbox:
            self.bbox = bbox

    def rotate90ccw(self):
        for c,l in self.layers:
            for i in range(len(l)):
                l[i] = ( -l[i][1], l[i][0], l[i][2] )
        self.bbox = rotatebbox90ccw( *self.bbox )

    def rotate90cw(self):
        for c,l in self.layers:
            for i in range(len(l)):
                l[i] = ( l[i][1], -l[i][0], l[i][2] )
        self.bbox = rotatebbox90cw( *self.bbox )

    def move(self,xoff,yoff):
        for c,l in self.layers:
            for i in range(len(l)):
                l[i] = ( l[i][0] + xoff, l[i][1] + yoff, l[i][2] )
        self.bbox = [
            self.bbox[0] + xoff,
            self.bbox[1] + yoff,
            self.bbox[2] + xoff,
            self.bbox[3] + yoff
        ]

    def in_bbox(self,x,y):
        if x < self.bbox[0] or x > self.bbox[2]:
            return False
        if y < self.bbox[1] or y > self.bbox[3]:
            return False
        return True

    def clip_to_bbox(self, remove = False):
        i = 0
        r = 0
        while i < len(self.layers):
            c,l = self.layers[i]
            j = 0
            a = True
            while j < len(l):
                x,y,d = l[j]
                if not self.in_bbox(x,y):
                    if remove:
                        del l[j]
                        r += 1
                    else:
                        x = min(max( self.bbox[0], x ), self.bbox[2] )
                        y = min(max( self.bbox[1], y ), self.bbox[3] )
                        l[j] = (x,y,d)
                        j += 1
                else:
                    a = False
                    j += 1
            if len(l) == 0 or a:
                r += len(l)
                del self.layers[i]
            else:
                i += 1
        return r

    def x_MO(self,unit):
        self.units = unit.lower()

    def x_FS(self,_format):
        m = re.match(r"(L|T)(A|I)X(\d\d)Y(\d\d)", _format)
        if m.group(1) == "T":
            self.z_tail_trim = True
        elif m.group(2) == "I":
            raise Exception("Incremental mode not implemented")
        self.f_x = tuple([ int(c) for c in m.group(3) ])
        self.f_y = tuple([ int(c) for c in m.group(4) ])

    def x_AM(self, name, desc):
        self.ap_macros[name] = desc

    def x_AD(self, desc):
        m = re.match(r"^D(\d+)(.+)$", desc)
        self.ap_descs[ int(m.group(1)) ] = m.group(2)

    def g_G36(self):
        self.layers.append( ("G36",[]) )

    def g_G37(self):
        self.layers.append( ("G37",[]) )

    def m_M02(self):
        pass

class Excellon():
    def __init__(self, fname):
        self.fname = fname
        self.tools = {}
        self.holes = []
        self.bbox = None

        self.units = None
        self.z_tail = True
        self.fmt = None

        self.readinghead = True
        self.read()

        if self.units and self.fmt:
            if len(self.tools) + len(self.holes) > 0:
                print("%s: OK" % fname)
            else:
                print("%s: empty" % fname)
        else:
            print(self)
            raise Exception("Parameters missing")

    def __str__(self):
        return "{1}:  Units: {0.units}, Tailing zeroes: {0.z_tail}, Format: {0.fmt}, Tools: {2}, Holes: {3}".format(
            self, self.fname, len(self.tools), sum( [ len(h[1]) for h in self.holes ] ) )

    def write(self,f):
        pass

    def read(self):
        f = open(self.fname,"rt")
        l = f.readline()
        while l:
            l = l.strip().upper()
            if l.startswith(";"):
                pass
            elif self.readinghead:
                if l == "M48":
                    pass
                elif l == "%" or l == "M95":
                    self.readinghead = False
                elif l.count(",") > 0:
                    l = l.split(",")
                    if hasattr(self,"h_%s" % l[0]):
                        getattr(self,"h_%s" % l[0])( *l[1:] )
                    else:
                        raise Exception("%s unknown: %s" % (l[0], ",".join(l)))
                elif l.startswith("T"):
                    m = re.match(r"^T(\d+)(.+)$", l)
                    self.tools[ int(m.group(1)) ] = m.group(2)
                else:
                    raise Exception("Unknown header command: %s" % l)
            elif l.startswith("G"):
                if hasattr(self,"g_%s" % l):
                    getattr(self,"g_%s" % l)()
                else:
                    raise Exception("G code %s unknown." % l)
            elif l.startswith("M"):
                if hasattr(self,"m_%s" % l):
                    getattr(self,"m_%s" % l)()
                else:
                    raise Exception("M code %s unknown." % l)
            elif l.startswith("T"):
                if int(l[1:]) in self.tools:
                    self.holes.append( (l,[]) )
            elif l.startswith("X"):
                self.add_coord( l )
            else:
                raise Exception("Unknown command: %s" % l)
            l = f.readline()
        f.close()

    def add_coord(self, coord):
        for c in coord:
            if c not in "XYD0123456789-":
                raise Exception("Unknown coordinate:", coord)
        x,_,y = coord.partition("Y")
        x = x.strip("X")
        if self.z_tail:
            x = int(x) / (10**self.fmt[1])
            y = int(y) / (10**self.fmt[1])
        else:
            Exception("Format without tailing zeroes not implemented")
        self.holes[-1][1].append( ( x,y ) )

    def calc_bbox(self):
        bbox = [None,None,None,None]
        for tool,holes in self.holes:
            for x,y in holes:
                if bbox[0] is None or x < bbox[0]:
                    bbox[0] = x
                if bbox[2] is None or x > bbox[2]:
                    bbox[2] = x
                if bbox[1] is None or y < bbox[1]:
                    bbox[1] = y
                if bbox[3] is None or y > bbox[3]:
                    bbox[3] = y
        if None not in bbox:
            self.bbox = bbox

    def rotate90ccw(self):
        for tool,holes in self.holes:
            for i in range(len(holes)):
                holes[i] = ( -holes[i][1], holes[i][0] )
        self.bbox = rotatebbox90ccw( *self.bbox )

    def rotate90cw(self):
        for tool,holes in self.holes:
            for i in range(len(holes)):
                holes[i] = ( holes[i][1], -holes[i][0] )
        self.bbox = rotatebbox90cw( *self.bbox )

    def move(self,xoff,yoff):
        for tool,holes in self.holes:
            for i in range(len(holes)):
                holes[i] = ( holes[i][0] + xoff, holes[i][1] + yoff )
        self.bbox = [
            self.bbox[0] + xoff,
            self.bbox[1] + yoff,
            self.bbox[2] + xoff,
            self.bbox[3] + yoff
        ]

    def h_FMAT(self,ver):
        if ver != "2":
            raise Exception("Unsupported format version: %s" % ver)

    def h_ICI(self,onoff):
        if onoff != "OFF":
            raise Exception("Incremental mode not supported.")

    def h_METRIC(self, zeroes, fmt):
        self.units = "mm"
        self.z_tail = zeroes == "TZ"
        before,_,after = fmt.partition(".")
        self.fmt = ( len(before), len(after) )

    def g_G90(self):
        pass

    def m_M71(self):
        if self.units != "mm":
            raise Exception("Mismatch between header units %s and M71" % self.units)

    def m_M30(self):
        pass

class Board():
    FILEEXTS = {
        "holes":[".XLN"],
        "layer-bottom":[".GBL"],
        "layer-top":[".GTL"],
        "mask-bottom":[".GBS"],
        "mask-top":[".GTS"],
        "outline":[".GKO"],
        "paste-bottom":[".GBP"],
        "paste-top":[".GTP"],
        "silk-bottom":[".GBO"],
        "silk-top":[".GTO"],
    }
    def __init__(self,filepattern):
        files = glob.glob(filepattern)
        self.files = {}
        for f in files:
            for l,exts in self.FILEEXTS.items():
                for ext in exts:
                    if f.endswith(ext):
                        if l in self.files:
                            print("Multiple files for Layer %s: %s, %s" % (l,self.files[l], f))
                        else:
                            self.files[l] = f
        self.layers = {}
        self.bbox = None
        for l,f in self.files.items():
            print("%s:" % l)
            if l in ["holes"]:
                self.layers[l] = Excellon(f)
            else:
                self.layers[l] = Gerber_RS_274X(f)
            self.layers[l].calc_bbox()
        if len(self.layers) > 0:
            self.calc_bbox()

    def save(self, dest):
        for n,l in self.layers.items():
            f = open( "%s%s" % (dest, self.FILEEXTS[n][0]), "w")
            l.write(f)
            f.close()

    def calc_bbox(self):
        if "outline" in self.layers and self.layers["outline"].bbox is not None:
            self.bbox = self.layers["outline"].bbox
        bbox = [None,None,None,None]
        wlayers = set()
        for n,l in self.layers.items():
            if l.bbox is not None:
                for i in [0,1]:
                    if bbox[i] is None or l.bbox[i] < bbox[i]:
                        bbox[i] = l.bbox[i]
                    if self.bbox is not None and l.bbox[i] < self.bbox[i]:
                        wlayers.add(n)
                for i in [2,3]:
                    if bbox[i] is None or l.bbox[i] > bbox[i]:
                        bbox[i] = l.bbox[i]
                    if self.bbox is not None and l.bbox[i] > self.bbox[i]:
                        wlayers.add(n)
        if self.bbox is None:
            bbox[0] -= 5
            bbox[1] -= 5
            bbox[2] += 5
            bbox[3] += 5
            self.bbox = bbox
            print("No outline layer, padding bbox by 5")
        elif self.bbox != bbox:
            print("Layers %s have Elements outside of outline" % ( ", ".join([l for l in wlayers]) ))
        for n,l in self.layers.items():
            l.bbox = self.bbox
        for l in wlayers:
            r = self.layers[l].clip_to_bbox()
            print("Clipped layer %s, removed %d coordinates." % ( l,r ))
        self.move( -self.bbox[0], -self.bbox[1] )

    def move(self, xoff, yoff):
        for n,l in self.layers.items():
            l.move( xoff, yoff )
        self.bbox = [
            self.bbox[0] + xoff,
            self.bbox[1] + yoff,
            self.bbox[2] + xoff,
            self.bbox[3] + yoff
        ]

    def rotate90ccw(self):
        for n,l in self.layers.items():
            l.rotate90ccw()
        self.bbox = rotatebbox90ccw( *self.bbox )
        self.move( -self.bbox[0], -self.bbox[1] )

    def rotate90cw(self):
        for n,l in self.layers.items():
            l.rotate90cw()
        self.bbox = rotatebbox90cw( *self.bbox )
        self.move( -self.bbox[0], -self.bbox[1] )

def rotatebbox90ccw( x1, y1, x2, y2 ):
    return [ -y2, x1, -y1, x2 ]

def rotatebbox90cw( x1, y1, x2, y2 ):
    return [ y1, -x2, y2, -x1 ]

for d in sys.argv[2:]:
    b = Board( os.path.join(d,"*") )
    b.rotate90ccw()
    b.save(sys.argv[1])
