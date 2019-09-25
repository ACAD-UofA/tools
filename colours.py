
#http://alumni.media.mit.edu/~wad/color/numbers.html
_colours="""
Black
    Lab: 0, 0, 0
    HSB: 0, 0, 0
    RGB: 0, 0, 0 
Dk. Gray
    Lab: 50, 0, 0
    HSB: 0, 0, 34
    RGB: 87, 87, 87 
Red
    Lab: 50, 58, 40
    HSB: 0, 80, 68
    RGB: 173, 35, 35 
Blue
    Lab: 50, 40, -77
    HSB: 229, 80, 95
    RGB: 42, 75, 215 
Green
    Lab: 50, -50, 42
    HSB: 114, 80, 41
    RGB: 29, 105, 20 
Brown
    Lab: 50, 20, 44
    HSB: 28, 80, 51
    RGB: 129, 74, 25 
Purple
    Lab: 50, 65, -61
    HSB: 275, 80, 75
    RGB: 129, 38, 192 
Lt. Gray
    Lab: 75, 0, 0
    HSB: 0, 0, 63
    RGB: 160, 160, 160 
Lt. Green
    Lab: 80, -34, 25
    HSB: 114, 38, 77
    RGB: 129, 197, 122 
Lt. Blue
    Lab: 80, 9, -34
    HSB: 229, 38, 100
    RGB: 157, 175, 255 
Cyan
    Lab: 80, -43, -14
    HSB: 180, 80, 82
    RGB: 41, 208, 208 
Orange
    Lab: 80, 28, 62
    HSB: 28, 80, 100
    RGB: 255, 146, 51 
Yellow
    Lab: 96, -19, 77
    HSB: 60, 80, 97
    RGB: 255, 238, 51 
Tan
    Lab: 91, 5, 12
    HSB: 28, 20, 93
    RGB: 233, 222, 187 
Pink
    Lab: 91, 15, 6
    HSB: 0, 20, 100
    RGB: 255, 205, 243 
White
    Lab: 100, 0, 0
    HSB: 0, 0, 100
    RGB: 255, 255, 255
"""

def _pal16():
    col = None
    for line in _colours.splitlines():
        if not line:
            continue
        if line[0] != " ":
            col = line
            continue
        line = line.strip()
        if line.startswith("RGB"):
            fields = map(int, line[5:].split(", "))
            rgb = "#{:02x}{:02x}{:02x}".format(*fields)
            yield col, rgb

pal16 = list(a for _,a in _pal16())

if __name__ == "__main__":
    print(pal16)
