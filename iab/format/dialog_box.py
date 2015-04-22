import markdown2
from IPython.display import HTML


class MessageBox(object):

    def __init__(self):
        self.boxes = []

    def add_warning(self, warning_note):
        box = """
        <div class='messageBox warning'>

        <header>
        <h1><i class="fa fa-exclamation-triangle"></i>
        <span style="clear:both">Warning!</span></h1>
        </header>

        <section>
        %s
        </section>


        </div>
        """

        self.boxes.append(box % markdown2.markdown(warning_note))
        return self

    def add_link(self, link_note):
        box = """
        <div class='messageBox link'>

        <header>
        <h1><i class="fa fa-link"></i>
        <span style="clear:both">Link for Additional Resources</span></h1>
        </header>

        <section>
        %s
        </section>


        </div>
        """

        self.boxes.append(box % markdown2.markdown(link_note))
        return self

    def add_developer(self, developer_note):
        box = """
        <div class='messageBox developer'>

        <header>
        <h1><i style="font-weight:bold" class="fa fa-terminal"></i>
        <span style="clear:both">Developer Note</span></h1>
        </header>

        <section>
        %s
        </section>


        </div>
        """

        self.boxes.append(box % markdown2.markdown(developer_note))
        return self

    def add_additional_info(self, general_note):
        box = """
        <div class='messageBox general'>

        <header>
        <h1><i class="fa fa-info-circle"></i></i>
        <span style="clear:both">Additional Information</span></h1>
        </header>

        <section>
        %s
        </section>


        </div>
        """

        self.boxes.append(box % markdown2.markdown(general_note))
        return self

    def _repr_html_(self):

        prefix = """
        <head>
        <link rel="stylesheet"href="//maxcdn.bootstrapcdn.com
        /font-awesome/4.3.0/css/font-awesome.min.css">
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

        .messageBox.warning header {
            background-color:#FFCC00;
            color:darkred;
        }

        .messageBox.warning section {
            background:#FFF9E5;
            color:darkred;
        }

        .messageBox.link header {
            background-color:dodgerblue;
            color:#fff;
        }

        .messageBox.link section {
            background:#e8f3ff;
            color:dodgerblue;
        }

        .messageBox.developer header {
            background-color:#000;
            color:#76EE00;
        }

        .messageBox.developer section {
            background:#e5e5e5;
            color:#000;
        }

        .messageBox.general header {
            background-color:#590059;
            color:#fff;
        }

        .messageBox.general section {
            background:#eee5ee;
            color:#590059;
        }


        </style>
        </head>

        """

        postfix = """

        """

        return "".join([prefix] + self.boxes + [postfix])


def make_box(section_text='', header_background='',
             header_text='', header_text_color='',
             icon='', section_background_color='',
             section_text_color='', style=''):
    return HTML("""
    <head>
    <link rel="stylesheet" href=
    "//maxcdn.bootstrapcdn.com/font-awesome/4.3.0/css/font-awesome.min.css">
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

    .messageBox.%s header {
        background-color:%s;
        color:%s;
    }

    .messageBox.%s section {
        background:%s;
        color:%s;
    }


    </style>
    </head>

    <div class='messageBox %s'>

    <header>
    <h1>%s    <span style="clear:both">%s</span></h1>
    </header>

    <section>
    %s
    </section>


    </div>
    """ % (style, header_background, header_text_color,
           style, section_background_color, section_text_color,
           style, icon, header_text, markdown2.markdown(section_text)))


def link(section_text,
         header_background='dodgerblue',
         header_text='Additional Resources',
         header_text_color='#fff',
         icon='<i class="fa fa-link"></i>',
         section_background_color='#e8f3ff',
         section_text_color='dodgerblue'):
    style = 'link'
    return make_box(section_text, header_background,
                    header_text, header_text_color,
                    icon, section_background_color,
                    section_text_color, style)


def warning(section_text,
            header_background='#FFCC00',
            header_text='Warning!',
            header_text_color='darkred',
            icon='<i class="fa fa-exclamation-triangle"></i>',
            section_background_color='#FFF9E5',
            section_text_color='darkred'):
    style = "warning"
    return make_box(section_text, header_background,
                    header_text, header_text_color,
                    icon, section_background_color,
                    section_text_color, style)


def additional_info(section_text,
                    header_background='#590059',
                    header_text='Additional Information',
                    header_text_color='#fff',
                    icon='<i class="fa fa-info-circle"></i></i>',
                    section_background_color='#eee5ee',
                    section_text_color='#590059'):
    style = 'additional_info'
    return make_box(section_text, header_background,
                    header_text, header_text_color,
                    icon, section_background_color,
                    section_text_color, style)


def developer_note(section_text,
                   header_background='#000',
                   header_text='Developer Note',
                   header_text_color='#76EE00',
                   icon="""<i style="font-weight:bold"
                   class="fa fa-terminal"></i>""",
                   section_background_color='#e5e5e5',
                   section_text_color='#000'):
    style = 'developer_note'
    return make_box(section_text, header_background,
                    header_text, header_text_color,
                    icon, section_background_color,
                    section_text_color, style)
