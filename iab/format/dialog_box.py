import markdown2
from IPython.display import HTML


def make_box(section_text='', header_background_color='none',
             header_text='', header_text_color='none',
             icon='', section_background_color='none',
             section_text_color='none', style=''):
    """Generic fucntion to create dialog box.

    This is a generic function that displays a string as HTML that is
    displayed inline in an IPython notebook. The HTML is automatically
    displayed when the function is called.

    Parameters
    ----------
    section_text : str
        The text that will be displayed in the main section of the dialog
        box
    header_background_color : str
        The color of the header background.
    header_text : str
        The text to be displayed in the header.
    header_text_color : str
        The color of the text in the header.
    icon : str
        The tag for the image or icon that will be displayed next to the
        header text.
    section_background_color : str
        The background color of the main text section.
    section_text_color : str
        The color of the main section text.
    style : str
        The name of the style used in the css class. This prevents the dialog
        box from overwriting other HTML objects and should therefore be
        unique.

    Returns
    -------
    IPython.core.display.HTML
        IPython object that displays the HTML inline in the IPython notebook

    Examples
    --------
    >>> from iab.format.dialog_box import make_box
    >>> make_box('foo',\
                 header_background_color='#000',\
                 header_text='Developer Note',\
                 header_text_color='#76EE00',\
                 icon="<i style='font-weight:bold'class\
                 ='fa fa-terminal'></i>",\
                 section_background_color='#e5e5e5',\
                 section_text_color='#000',\
                 style='developer_note')\
    <IPython.core.display.HTML object>

    .. shownumpydoc

    """
    return HTML("""
    <head>
    <link rel="stylesheet" href="//maxcdn.bootstrapcdn.com/
    font-awesome/4.3.0/css/font-awesome.min.css">
    <style>

    .messageBox h1 {
        margin:3px 0px 5px 0px;
    }

    .messageBox header {
       border-top-left-radius:5px;
       border-top-right-radius:5px;
        padding-left:15px;
        padding-top:1px;
        padding-bottom:1px;
        line-height:10px;

    }

    .messageBox section {
        padding:15px;

    }

    .messageBox.%(style)s header {
        background-color:%(header_background_color)s;
        color:%(header_text_color)s;
    }

    .messageBox.%(style)s section {
        background:%(section_background_color)s;
        color:%(section_text_color)s;
    }


    </style>
    </head>

    <div class='messageBox %(style)s'>

    <header>
    <h1>%(icon)s    <span style="clear:both">%(header_text)s</span></h1>
    </header>

    <section>
    %(section_text)s
    </section>

    </div>
    """ % {"style": style,
           "header_background_color": header_background_color,
           "header_text_color": header_text_color,
           "section_background_color": section_background_color,
           "section_text_color": section_text_color,
           "icon": icon,
           "header_text": header_text,
           "section_text": markdown2.markdown(section_text)})


def link(section_text):
    return make_box(section_text=section_text,
                    header_background_color='dodgerblue',
                    header_text='Additional Resources',
                    header_text_color='#fff',
                    icon='<i class="fa fa-link"></i>',
                    section_background_color='#e8f3ff',
                    section_text_color='dodgerblue',
                    style='link_box')


def warning(section_text):
    return make_box(section_text=section_text,
                    header_background_color='#FFCC00',
                    header_text='Warning!',
                    header_text_color='darkred',
                    icon='<i class="fa fa-exclamation-triangle"></i>',
                    section_background_color='#FFF9E5',
                    section_text_color='darkred',
                    style="warning_box")


def additional_info(section_text):
    return make_box(section_text=section_text,
                    header_background_color='#590059',
                    header_text='Additional Information',
                    header_text_color='#fff',
                    icon='<i class="fa fa-info-circle"></i></i>',
                    section_background_color='#eee5ee',
                    section_text_color='#590059',
                    style='additional_box')


def developer_note(section_text):
    return make_box(section_text=section_text,
                    header_background_color='#000',
                    header_text='Developer Note',
                    header_text_color='#76EE00',
                    icon="""<i style="font-weight:bold"
                    class="fa fa-terminal"></i>""",
                    section_background_color='#e5e5e5',
                    section_text_color='#000',
                    style='developer_note')
