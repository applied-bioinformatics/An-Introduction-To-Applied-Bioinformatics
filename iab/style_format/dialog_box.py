import markdown2

class MessageBox(object):

    def __init__(self):
        self.boxes = []

    def add_warning(self, warning_note):
        box = """
        <div class='messageBox warning'>

        <header>
        <h1><i class="fa fa-exclamation-triangle"></i>    <span style="clear:both">Warning!</span></h1>
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
        <h1><i class="fa fa-link"></i>    <span style="clear:both">Link for Additional Resources</span></h1>
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
        <h1><i style="font-weight:bold" class="fa fa-terminal"></i>    <span style="clear:both">Developer Note</span></h1>
        </header>

        <section>
        %s
        </section>


        </div>
        """

        self.boxes.append(box % markdown2.markdown(developer_note))
        return self

    def add_general(self, general_note):
        box = """
        <div class='messageBox general'>

        <header>
        <h1><i class="fa fa-info-circle"></i></i>    <span style="clear:both">Additional Information</span></h1>
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
        <link rel="stylesheet" href="//maxcdn.bootstrapcdn.com/font-awesome/4.3.0/css/font-awesome.min.css">
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
