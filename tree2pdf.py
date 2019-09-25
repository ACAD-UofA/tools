#!/usr/bin/env python

from __future__ import print_function
import sys
import ete3
from PyQt4 import QtCore
from PyQt4.QtGui import QGraphicsTextItem, QGraphicsRectItem, QPen, QFont


def name2date(name):
    last = name.split("_")[-1]
    try:
        year = int(last)
    except ValueError:
        year = 0
    return year

def isACAD(name):
    if name[0] == "A" or "SW18" in name or "BB20087" in name or "Q229" in name:
        return True
    else:
        return False

def mapname(name):

    if "camelopardalis" in name:
        return ("Giraffa camelopardalis", None, "giraffe")
    if "bubalis" in name:
        return ("Bubalus bubalis", None, "water buffalo")
    if "grunniens" in name:
        return ("Bos grunniens", None, "yak")
    if "frontalis" in name:
        return ("Bos frontalis", None, "gayal")
    if "javanicus" in name:
        return ("Bos javanicus", None, "banteng")
    if "caucasicus" in name or "bonasus" in name:
        return ("Bison bonasus", None, "European bison")
#    if "caucasicus" in name:
#        return ("Bison bonasus", "caucasicus", "European bison (highland)")
#    if "bonasus" in name:
#        return ("Bison bonasus", "bonasus", "European bison (lowland)")
    if "priscus" in name:
        return ("Bison priscus", None, "steppe bison")
    if "bison" in name:
        return ("Bison bison", None, "American bison")
    if "taurus" in name:
        if "Gir" in name or "Nelore" in name or "Brahman" in name:
            return ("Bos taurus", "indicus", "zebu")
        else:
            return ("Bos taurus", "taurus", "cattle")
    if "primigenius" in name:
        return ("Bos primigenius", None, "aurochs")
    if "aries" in name:
        return ("Ovis aries", None, "sheep")

    # unknown
    return (name, name)

def trinomial_face(binom, spp, ftype="Verdana", fsize=10,
                 fgcolor="black", fstyle="normal",
                 bold=False):
    """
    Stolen from:
    https://connorskennerton.wordpress.com/2016/11/16/python-ete3-formatting-organism-names-the-way-i-want/
    """

    font = QFont(ftype, fsize)
    font.setBold(bold)
    if fstyle == "italic":
        font.setStyle(QFont.StyleItalic)
    elif fstyle == "oblique":
        font.setStyle(QFont.StyleOblique)

    text = QGraphicsTextItem()
    text.setFont(font)
    if spp is None:
        text.setHtml("<i>{}</i>".format(binom))
    else:
        text.setHtml("<i>{}</i> {}".format(binom, spp))

    rect = QGraphicsRectItem(0,0,text.boundingRect().width(),text.boundingRect().height()-10)
    text.setParentItem(rect)
    center = rect.boundingRect().center()
    text.setPos(rect.boundingRect().x(),
        center.y() - text.boundingRect().height()/2)
 
    # no border
    rect.setPen(QPen(QtCore.Qt.NoPen))

    return ete3.faces.StaticItemFace(rect)

def bovid_layout(node, fgcolor, labelnodes):

    node.img_style["hz_line_width"] = node.img_style["vt_line_width"] = 3
    mar = 5

    if node.is_leaf():
        if isACAD(node.name):
            node.img_style["size"] = 6
            node.img_style["shape"] = "square"
            node.img_style["fgcolor"] = "red"
        else:
            node.img_style["size"] = 0
            node.img_style["shape"] = "square"

        n_face = ete3.faces.TextFace(node.name, fsize=10)
        n_face.margin_left = n_face.margin_right = mar
        ete3.faces.add_face_to_node(n_face, node, 0, position="branch-right")
    else:
        node.img_style["size"] = 0
        node.img_style["shape"] = "square"
        node.img_style["fgcolor"] = fgcolor

    if node in labelnodes:
        binom, spp, common = labelnodes[node]

        common_face = ete3.faces.TextFace(common, fsize=16, bold=True, fgcolor=fgcolor)
        common_face.margin_right = common_face.margin_left = mar
        ete3.faces.add_face_to_node(common_face, node, column=0, position="aligned")

        trinom_face = trinomial_face(binom, spp, fsize=10)
        ete3.faces.add_face_to_node(trinom_face, node, column=0, position="aligned")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("usage: {} in.nwk out.pdf".format(sys.argv[0]), file=sys.stderr)
        exit(1)

    nwk_fn = sys.argv[1]
    pdf_fn = sys.argv[2]

    tree = ete3.Tree(nwk_fn)

    # midpoint rooting
    midpoint = tree.get_midpoint_outgroup()
    tree.set_outgroup(midpoint)

    # order children nodes by the number of decendent nodes
    for node in tree.traverse():
        if len(node.children) != 2:
            continue
        if len(node.children[0].get_leaves()) > len(node.children[1].get_leaves()):
            node.swap_children()

    for node in tree.iter_leaves():
        if node.name.endswith("_0"):
            node.name = node.name[:-2]

    pal16 = ['#000000', '#575757', '#ad2323', '#2a4bd7',
            '#1d6914', '#814a19', '#8126c0', '#a0a0a0',
            '#81c57a', '#9dafff', '#29d0d0', '#ff9233',
            '#ffee33', '#e9debb', '#ffcdf3', '#ffffff']
    pal2 = ['#a0a0a0', '#e9debb']*4+['#a0a0a0']+['#a0a0a0', '#e9debb']*4
    pal2 = ['#a0a0a0', '#e9debb']*8
    bgcol = iter(pal16[2:])
    bgcol = iter(pal2)

    nodestyles = {}
    for node in tree.iter_leaves():
        binom, spp, common = mapname(node.name)
        if (binom,spp) in nodestyles:
            nodestyles[(binom,spp)][1].append(node)
            continue
        col = next(bgcol)
        s = ete3.NodeStyle()
        s["bgcolor"] = col
        #s["hz_line_color"] = s["vt_line_color"] = col
        #s["hz_line_width"] = s["vt_line_width"] = 1
        nodestyles[(binom,spp)] = (s, [node])

    import copy

    labelnodes = {}
    for (binom,spp), (ns, nodelist) in nodestyles.iteritems():
        if len(nodelist) == 1:
            anc = nodelist[0]
        else:
            anc = tree.get_common_ancestor(*nodelist)
        #labelnodes[nodelist[len(nodelist)/2]] = mapname(nodelist[0].name)
        labelnodes[nodelist[0]] = mapname(nodelist[0].name)

#        for n in anc.iter_descendants():
        for n in nodelist:
            ns2 = copy.deepcopy(ns)
            n.set_style(ns2)
        anc.set_style(ns)

    fgcol = "black"

    style = ete3.TreeStyle()
    style.mode = "r"
    style.layout_fn = lambda n: bovid_layout(n, fgcol, labelnodes)
    style.show_leaf_name = False
    mar = 10
    style.margin_right = style.margin_left = mar
    style.margin_top = mar
    #style.margin_bottom = mar
    style.title.add_face(ete3.faces.TextFace("Nuclear phylogeny (pairwise distance with neighbour-joining)", fsize=20), column=0)
    style.scale = 3*1920
#    style.branch_vertical_margin = 10

    #tree.render(pdf_fn, w=1920, h=1080, units="px", tree_style=style)
    tree.render(pdf_fn, h=1080, tree_style=style)
