from markdown.blockprocessors import BlockProcessor
from MooseCommonExtension import MooseCommonExtension
import re
from markdown.util import etree

class MooseCSS(BlockProcessor, MooseCommonExtension):
  """
  Markdown extension for applying CSS styles to paragraph
  Markdown syntax is:
   !css <options>
   Paragraph text here

  Where <options> are key=value pairs.
  """

  RE = re.compile(r'^!\ ?css(.*)')
  # If there are multiple css blocks on the same page then
  # they need to have different ids
  MATCHES_FOUND = 0

  def __init__(self, parser, **kwargs):
   MooseCommonExtension.__init__(self, **kwargs)
   BlockProcessor.__init__(self, parser)

  def test(self, parent, block):
   """
   Test to see if we should process this block of markdown.
   Inherited from BlockProcessor.
   """
   return self.RE.search(block)

  def run(self, parent, blocks):
   """
   Called when it is determined that we can process this block.
   This will convert the markdown into HTML
   """
   sibling = self.lastChild(parent)
   block = blocks.pop(0)
   m = self.RE.search(block)

   if m:
    # Parse out the options on the css line
    options, styles = self.getSettings(m.group(1))
    block = block[m.end() + 1:] # removes the css line

   block, paragraph = self.detab(block)
   if m:
    top_div = etree.SubElement(parent, 'div')
    self.createCSS(top_div, styles, paragraph)
   else:
    top_div = sibling

   self.parser.parseChunk(top_div, block)

  def createCSS(self, top_div, styles, paragraph):
    """
    Creates the actual HTML required for the CSS paragraph to work.
    Input:
     top_div: div element that will contain the paragraph element
     styles[dict]: The CSS style attributes
     paragraph: the actual text within the <p></p> element
    """

    p_el = self.addStyle(etree.SubElement(top_div, 'p'), **styles)
    p_el.text = paragraph
