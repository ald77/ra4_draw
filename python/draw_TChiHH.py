from pyfeyn.user import *
from math import sin,cos,tan
import pyx 

# paint styles
# https://pyfeyn.hepforge.org/doc/pyfeyn.paint-module.html

chi10=r"$\tilde{\chi}^0_1$"
goldstino=r"$\tilde{G}$"

susycolor = RED
# susycolor = BLACK
linewidth = THICK2

r_blob = 0.4

dx_1 = 1.5
ang_1 = 30./180.*3.1415

dx_h = 1.
ang_h = 40./180.*3.1415
dx_g = 1.7
ang_g = -5./180.*3.1415

processOptions()
fd = FeynDiagram()

in1  = Point(-2.2, -1.5)
in2  = Point(-2.2, 1.5)
blob = Circle(-1,0, radius=r_blob).setFillStyle(GRAY)
P1 = Fermion(in1, blob).addLabel("$p$",displace=-0.2)
P1_label = MultiLine(in1, blob,n=2, dist=0.1)
P2 = Fermion(in2, blob).addLabel("$p$",displace=.17)
P2_label = MultiLine(in2, blob,n=2, dist=0.1)

vtx_top =  Vertex(dx_1,      dx_1*tan(ang_1),  mark=CircleMark())
vtx_toph = Vertex(dx_1+dx_h, dx_1*tan(ang_1)+dx_h*tan(ang_h))
vtx_topg = Vertex(dx_1+dx_g, dx_1*tan(ang_1)+dx_g*tan(ang_g))

vtx_bot =  Vertex(dx_1,      -dx_1*tan(ang_1), mark=CircleMark())
vtx_both = Vertex(dx_1+dx_h, -(dx_1*tan(ang_1)+dx_h*tan(ang_h)))
vtx_botg = Vertex(dx_1+dx_g, -(dx_1*tan(ang_1)+dx_g*tan(ang_g)))


chi10_top = Vector(blob, vtx_top).setAmplitude(0.1).setFrequency(0.4)
chi10_topin = Line(blob, vtx_top).addStyle(linewidth).addStyle(susycolor)
# chi10_top = Gaugino(blob, vtx_top)
# chi10_top.set3D()
chi10_top.addLabel(chi10)
chi10_top.addStyle(linewidth).addStyle(susycolor)
h_top = Higgs(vtx_top, vtx_toph).addLabel("$h$",displace=0.05, pos = 1.3)
# h_top = Higgs(vtx_top, vtx_toph).addLabel("$h$",angle=90, pos = 1.3)
h_top.addStyle(linewidth)
g_top = Ghost(vtx_top, vtx_topg).addLabel(goldstino,displace=0.005, pos = 1.2)
g_top.addStyle(linewidth).addStyle(susycolor)

chi10_bot = Vector(blob, vtx_bot).setAmplitude(0.1).setFrequency(0.4)
chi10_botin = Line(blob, vtx_bot).addStyle(linewidth).addStyle(susycolor)
# chi10_bot = Gaugino(blob, vtx_bot)
# chi10_bot.set3D()
chi10_bot.addLabel(chi10,displace=.34)
chi10_bot.addStyle(linewidth).addStyle(susycolor)
chi10_bot.invert()

g_bot = Ghost(vtx_bot, vtx_botg).addLabel(goldstino,displace=0.01, pos = 1.2)
g_bot.addStyle(linewidth).addStyle(susycolor)
h_bot = Higgs(vtx_bot, vtx_both).addLabel("$h$",displace=0., pos = 1.3)
h_bot.addStyle(linewidth)


fd.draw("TChiHH_diag.pdf")

# os.system ("convert %s_feyn.pdf -transparent white %s_feyn.png" % (name, name))