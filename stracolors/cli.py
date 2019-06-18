import sys

import argh

from stracolors import commands as cmd

def main():
    parser = argh.ArghParser()
    parser.add_commands([cmd.call_palette,
                        cmd.call_palette_binary,
                        cmd.uoe_colors])
    parser.dispatch()

if __name__ == "__main__":
    sys.exit(main())
