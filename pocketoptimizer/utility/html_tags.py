"""Classes to generate HTML in Python


The HTMLTags module defines a class for all the valid HTML tags, written in
uppercase letters. To create a piece of HTML, the general syntax is :
    t = TAG(value, key1=val1,key2=val2,...)

so that "print t" results in :
    <tag key1="val1" key2="val2" ...>value</tag>

For instance :
    print A('bar', href="foo") ==> <A href="foo">bar</A>

To generate HTML attributes without value, give them the value True :
    print OPTION('foo',SELECTED=True,value=5) ==>
            <option value="5" SELECTED>

The value argument can be an instance of an HTML class, so that
you can nest tags, like this :
    print B(I('foo')) ==> <b><i>foo</i></b>
you can also provide several inner tags:
    print B(I('foo'),I('bar')) ==> <b><i>foo</i><i>bar</i></b>


For complex expressions, a tag can be nested in another using the method insert
Considering the HTML document as a tree, this means "add child" :

    form = FORM(action="foo")
    form.insert(INPUT(name="bar"))
    form.insert(INPUT(Type="submit",value="Ok"))

insert can add several childs at once:
    form = FORM(action="foo")
    form.insert(INPUT(name="bar"), INPUT(Type="submit",value="Ok"))


generates the rows of a table showing the squares of integers from 0 to 99

A simple document can be produced by :
    print HTML( HEAD(TITLE('Test document')),
        BODY(H1('This is a test document'),
             'First line',BR(),
             'Second line'))

This will be rendered as :
    <html>
     <head>
      <title>Test document</title>
     </head>
     <body>
      <h1>This is a test document</h1>
    First line
      <br>
    Second line
     </body>
    </html>

If the document is more complex it is more readable to create the elements
first, then to print the whole result in one instruction. For example :

head = HEAD()
head.insert(TITLE('Record collection'))
head.insert(LINK(rel="Stylesheet",href="doc.css"))

title = H1('My record collection')
table = TABLE()
table.insert(TR(TH('Title'), TH('Artist')))
for rec in records:
    row = TR()
    # note the attribute key Class with leading uppercase
    # because "class" is a Python keyword
    row.insert(TD(rec.title,Class="title"), TD(rec.artist,Class="artist"))
    table.insert(row)

print HTML(head, BODY(title, table))
"""

from io import StringIO
from typing import List, Union, NoReturn


class TAG:
    """Generic class for tags"""

    def __init__(self, *values: List[Union['TAG', str]], **attrs: List[Union[str, bool]]):
        """
        Constructor

        Parameters
        ----------
        values: list
            Values of the html tag. Can contain other tags or strings
        attrs: list
            Attributes of the html tag
        """
        self.tag = self.__class__.__name__
        self.attrs = attrs
        self.parent = None
        self.value = None
        self.children = []
        for value in values:
            if isinstance(value, TAG):
                self.insert(value)
            else:
                self.value = value

    def __str__(self) -> str:
        """
        Convert tag to string

        Returns
        -------
        String representation of tag
        """
        return self.to_string()

    def __mul__(self, n: int) -> 'TAG':
        """
        Replicate self n times, with tag first : TAG * n

        Parameters
        ----------
        n: int
            Number of replicate tags to be created

        Returns
        -------
        The multiplied tag
        """
        res = TAG()
        res.tag = self.tag
        res.value = self.value
        res.attrs = self.attrs
        for dummy in range(n - 1):
            res += self
        return res

    def __rmul__(self, n: int) -> 'TAG':
        """
        Replicate self n times, with n first : n * TAG

        Parameters
        ----------
        n: int
            Number of replicate tags to be created

        Returns
        -------
        The multiplied tag
        """
        return self * n

    def to_string(self, indent: int = 0) -> str:
        """
        Prints the html tag and its contents/subtags.
        Subtags are further indented.

        Parameters
        ----------
        indent: int
            Number of spaces that subtags are to be indented

        Returns
        -------
        String representing this tag
        """
        result_string = StringIO()
        w = result_string.write
        if self.tag != "TEXT":
            if self.parent is not None and self.parent.tag in ONE_LINE:
                w(' ' * indent)
            w("<%s" % self.tag)
            # attributes which will produce arg = "val"
            attr1 = [k for k in self.attrs
                     if not isinstance(self.attrs[k], bool)]
            w("".join([' %s="%s"'
                       % (k.replace('_', '-'), self.attrs[k]) for k in attr1]))
            # attributes with no argument
            # if value is False, don't generate anything
            attr2 = [k for k in self.attrs if self.attrs[k] is True]
            w("".join([' %s' % k for k in attr2]))
            if self.tag in NON_CLOSING_TAGS:
                w(" /")
            w(">")
        if self.tag in ONE_LINE:
            w('\n')
        if isinstance(self.value, str):
            if self.tag in ONE_LINE:
                w(' ' * indent)
            w(self.value)
        for child in self.children:
            w(child.to_string(indent + 1))
        if self.tag in CLOSING_TAGS:
            if self.tag in ONE_LINE:
                w(' ' * indent)
            w("</%s>" % self.tag)
        if self.tag in LINE_BREAK_AFTER:
            w('\n')
        return result_string.getvalue()

    def insert(self, *others: List[Union['TAG', str]]) -> NoReturn:
        """
        Add children to tag

        Parameters
        ----------
        others: list
            Children to be added (tags and/or value strings)
        """
        for other in others:
            if isinstance(other, str):
                other = TEXT(other)
            self.children.append(other)
            other.parent = self

    def add_attributes(self, **attrs: List[Union[str, bool]]) -> NoReturn:
        """
        Add attributes to this tag

        Parameters
        ----------
        attrs: list
            Attributes to be added
        """
        self.attrs.update(attrs)

    # list of tags, from the HTML 4.01 specification


CLOSING_TAGS = ['A', 'ABBR', 'ACRONYM', 'ADDRESS', 'APPLET',
                'B', 'BDO', 'BIG', 'BLOCKQUOTE', 'BUTTON',
                'CAPTION', 'CENTER', 'CITE', 'CODE',
                'DEL', 'DFN', 'DIR', 'DIV', 'DL',
                'EM', 'FIELDSET', 'FONT', 'FORM', 'FRAMESET',
                'H1', 'H2', 'H3', 'H4', 'H5', 'H6',
                'I', 'IFRAME', 'INS', 'KBD', 'LABEL', 'LEGEND',
                'MAP', 'MENU', 'NOFRAMES', 'NOSCRIPT', 'OBJECT',
                'OL', 'OPTGROUP', 'PRE', 'Q', 'S', 'SAMP',
                'SCRIPT', 'SELECT', 'SMALL', 'SPAN', 'STRIKE',
                'STRONG', 'STYLE', 'SUB', 'SUP', 'TABLE',
                'TEXTAREA', 'TITLE', 'TT', 'U', 'UL',
                'VAR', 'BODY', 'COLGROUP', 'DD', 'DT', 'HEAD',
                'HTML', 'LI', 'P', 'TBODY', 'OPTION',
                'TD', 'TFOOT', 'TH', 'THEAD', 'TR']

NON_CLOSING_TAGS = ['AREA', 'BASE', 'BASEFONT', 'BR', 'COL', 'FRAME',
                    'HR', 'IMG', 'INPUT', 'ISINDEX', 'LINK',
                    'META', 'PARAM']

# whitespace-insensitive tags, determines pretty-print rendering
LINE_BREAK_AFTER = NON_CLOSING_TAGS + ['HTML', 'HEAD', 'BODY',
                                       'FRAMESET', 'FRAME', 'TITLE', 'SCRIPT',
                                       'TABLE', 'TR', 'TD', 'COLGROUP', 'TH', 'SELECT', 'OPTION',
                                       'FORM', 'H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'P', 'THEAD', 'TBODY']

# tags whose opening tag should be alone in its line
ONE_LINE = ['HTML', 'HEAD', 'BODY', 'FRAMESET',
            'SCRIPT', 'TABLE', 'TR', 'SELECT', 'OPTION',
            'FORM', 'THEAD', 'TBODY']


class TEXT(TAG):
    pass


class A(TAG):
    pass


class ABBR(TAG):
    pass


class ACRONYM(TAG):
    pass


class ADDRESS(TAG):
    pass


class APPLET(TAG):
    pass


class B(TAG):
    pass


class BDO(TAG):
    pass


class BIG(TAG):
    pass


class BLOCKQUOTE(TAG):
    pass


class BUTTON(TAG):
    pass


class CAPTION(TAG):
    pass


class CENTER(TAG):
    pass


class CITE(TAG):
    pass


class CODE(TAG):
    pass


class DEL(TAG):
    pass


class DFN(TAG):
    pass


class DIR(TAG):
    pass


class DIV(TAG):
    pass


class DL(TAG):
    pass


class EM(TAG):
    pass


class FIELDSET(TAG):
    pass


class FONT(TAG):
    pass


class FORM(TAG):
    pass


class FRAMESET(TAG):
    pass


class H1(TAG):
    pass


class H2(TAG):
    pass


class H3(TAG):
    pass


class H4(TAG):
    pass


class H5(TAG):
    pass


class H6(TAG):
    pass


class I(TAG):
    pass


class IFRAME(TAG):
    pass


class INS(TAG):
    pass


class KBD(TAG):
    pass


class LABEL(TAG):
    pass


class LEGEND(TAG):
    pass


class MAP(TAG):
    pass


class MENU(TAG):
    pass


class NOFRAMES(TAG):
    pass


class NOSCRIPT(TAG):
    pass


class OBJECT(TAG):
    pass


class OL(TAG):
    pass


class OPTGROUP(TAG):
    pass


class PRE(TAG):
    pass


class Q(TAG):
    pass


class S(TAG):
    pass


class SAMP(TAG):
    pass


class SCRIPT(TAG):
    pass


class SELECT(TAG):
    pass


class SMALL(TAG):
    pass


class SPAN(TAG):
    pass


class STRIKE(TAG):
    pass


class STRONG(TAG):
    pass


class STYLE(TAG):
    pass


class SUB(TAG):
    pass


class SUP(TAG):
    pass


class TABLE(TAG):
    pass


class TEXTAREA(TAG):
    pass


class TITLE(TAG):
    pass


class TT(TAG):
    pass


class U(TAG):
    pass


class UL(TAG):
    pass


class VAR(TAG):
    pass


class BODY(TAG):
    pass


class COLGROUP(TAG):
    pass


class DD(TAG):
    pass


class DT(TAG):
    pass


class HEAD(TAG):
    pass


class HTML(TAG):
    pass


class LI(TAG):
    pass


class P(TAG):
    pass


class TBODY(TAG):
    pass


class OPTION(TAG):
    pass


class TD(TAG):
    pass


class TFOOT(TAG):
    pass


class TH(TAG):
    pass


class THEAD(TAG):
    pass


class TR(TAG):
    pass


class AREA(TAG):
    pass


class BASE(TAG):
    pass


class BASEFONT(TAG):
    pass


class BR(TAG):
    pass


class COL(TAG):
    pass


class FRAME(TAG):
    pass


class HR(TAG):
    pass


class IMG(TAG):
    pass


class INPUT(TAG):
    pass


class ISINDEX(TAG):
    pass


class LINK(TAG):
    pass


class META(TAG):
    pass


class PARAM(TAG):
    pass
