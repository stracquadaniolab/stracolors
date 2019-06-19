import sys

import argh
import matplotlib as mpl
mpl.use("TkAgg")

from matplotlib import pyplot as plt
plt.switch_backend("TkAgg")

from stracolors import commands as cmd

def main():
    parser = argh.ArghParser()
    parser.add_commands([cmd.call_palette,
                        cmd.call_palette_binary,
                        cmd.uoe_colors,
                        cmd.show_choices])
    parser.dispatch()

if __name__ == "__main__":
    sys.exit(main())
