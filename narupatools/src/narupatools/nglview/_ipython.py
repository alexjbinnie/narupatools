from IPython.core.magic import Magics, magics_class, line_magic

from ._show import show


@magics_class
class NGLMagics(Magics):

    @line_magic
    def ngl(self, line):
        """
        Display the given object using nglview.
        """
        return show(self.shell.ev(line))