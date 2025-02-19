from chimerax.core.tools import ToolInstance

class BasicPrivateerTool(ToolInstance):
    # Inheriting from ToolInstance makes us known to the ChimeraX tool mangager,
    # so we can be notified and take appropriate action when sessions are closed,
    # saved, or restored, and we will be listed among running tools and so on.
    #
    # If cleaning up is needed on finish, override the 'delete' method
    # but be sure to call 'delete' from the superclass at the end.

    SESSION_ENDURING = False    # Does this instance persist when session closes
    SESSION_SAVE = True         # We do save/restore in sessions
    help = "help:privateer_for_chimeraX_docs.html"
                                # Let ChimeraX know about our help page

    def __init__(self, session, tool_name):
        # 'session'   - chimerax.core.session.Session instance
        # 'tool_name' - string

        # Initialize base class.
        super().__init__(session, tool_name)

        # Set name displayed on title bar (defaults to tool_name)
        # Must be after the superclass init, which would override it.
        self.display_name = "Privateer - Validate Carbohydrates (save)"

        # Create the main window for our tool.  The window object will have
        # a 'ui_area' where we place the widgets composing our interface.
        # The window isn't shown until we call its 'manage' method.
        from chimerax.ui import MainToolWindow
        self.tool_window = MainToolWindow(self)

        # We will be adding an item to the tool's context menu, so override
        # the default MainToolWindow fill_context_menu method
        self.tool_window.fill_context_menu = self.fill_context_menu

        self._build_ui(session)

    def _build_ui(self,session):
        # Put our widgets in the tool window

        # We will use an editable single-line text input field (QLineEdit)
        # with a descriptive text label to the left of it (QLabel).  To
        # arrange them horizontally side by side we use QHBoxLayout
        from Qt.QtWidgets import QFormLayout, QComboBox, QPushButton
        from chimerax import atomic
        models = atomic.all_structures(session)
        layout = QFormLayout()
        label1 = "Model ID:"
        label2 = "Output Folder Path:"
        self.combo_box = QComboBox()
        for i, m in enumerate(models):
            self.combo_box.addItem(str(i+1))
        #self.line_edit = QLineEdit()
        self.run_button = QPushButton("Run Privateer")
        self.file_button = QPushButton("Select Output Directory")

        layout.addRow(label1,self.combo_box)
        layout.addRow(label2,self.file_button)
        layout.addWidget(self.run_button)

        layout.setFieldGrowthPolicy(QFormLayout.ExpandingFieldsGrow)

        self.file_button.clicked.connect(self.selectDirectoryDialog)

        # Arrange for our 'return_pressed' method to be called when the
        # user presses the Return key
        self.run_button.clicked.connect(self.button_pressed)

        # Set the layout as the contents of our window
        self.tool_window.ui_area.setLayout(layout)

        # Show the window on the user-preferred side of the ChimeraX
        # main window
        self.tool_window.manage('side')
    
    def selectDirectoryDialog(self):
        from Qt.QtWidgets import QFileDialog
        self.file_dialog = QFileDialog()
        self.file_dialog.setWindowTitle("Select Output Directory")
        self.file_dialog.setFileMode(QFileDialog.FileMode.Directory)
        if self.file_dialog.exec():
            selected_directory = self.file_dialog.selectedFiles()[0]
            self.file_button.setText(str(selected_directory))

    def button_pressed(self):
        # The user has pressed the Return key; run the privateer command using their inputs
        from chimerax.core.commands import run
        # ToolInstance has a 'session' attribute...
        run(self.session, f"privateer_validation {self.file_dialog.selectedFiles()[0]} {self.combo_box.currentText()}") 

    def fill_context_menu(self, menu, x, y):
        # Add any tool-specific items to the given context menu (a QMenu instance).
        # The menu will then be automatically filled out with generic tool-related actions
        # (e.g. Hide Tool, Help, Dockable Tool, etc.) 
        #
        # The x,y args are the x() and y() values of QContextMenuEvent, in the rare case
        # where the items put in the menu depends on where in the tool interface the menu
        # was raised.
        from Qt.QtGui import QAction
        clear_action = QAction("Clear", menu)
        clear_action.triggered.connect(lambda *args: self.file_dialog.clear())
        clear_action.triggered.connect(lambda *args: self.file_button.clear())
        clear_action.triggered.connect(lambda *args: self.combo_box.clear())
        clear_action.triggered.connect(lambda *args: self.run_button.clear())
        menu.addAction(clear_action)

class FancyPrivateerTool(ToolInstance):
    # Inheriting from ToolInstance makes us known to the ChimeraX tool mangager,
    # so we can be notified and take appropriate action when sessions are closed,
    # saved, or restored, and we will be listed among running tools and so on.
    #
    # If cleaning up is needed on finish, override the 'delete' method
    # but be sure to call 'delete' from the superclass at the end.

    SESSION_ENDURING = False    # Does this instance persist when session closes
    SESSION_SAVE = True         # We do save/restore in sessions
    help = "help:privateer_for_chimeraX_docs.html"
                                # Let ChimeraX know about our help page

    def __init__(self, session, tool_name):
        # 'session'   - chimerax.core.session.Session instance
        # 'tool_name' - string

        # Initialize base class.
        super().__init__(session, tool_name)

        # Set name displayed on title bar (defaults to tool_name)
        # Must be after the superclass init, which would override it.
        self.display_name = "Privateer - Validate Carbohydrates (display)"

        # Create the main window for our tool.  The window object will have
        # a 'ui_area' where we place the widgets composing our interface.
        # The window isn't shown until we call its 'manage' method.
        from chimerax.ui import MainToolWindow
        self.tool_window = MainToolWindow(self)

        # We will be adding an item to the tool's context menu, so override
        # the default MainToolWindow fill_context_menu method
        self.tool_window.fill_context_menu = self.fill_context_menu

        self._build_ui()

    def _build_ui(self):
        # Put our widgets in the tool window

        # We will use an editable single-line text input field (QLineEdit)
        # with a descriptive text label to the left of it (QLabel).  To
        # arrange them horizontally side by side we use QHBoxLayout
        from Qt.QtWidgets import QFormLayout, QComboBox, QPushButton
        from chimerax import atomic
        models = atomic.all_structures(self.session)
        layout = QFormLayout()
        label1 = "Model ID:"
        self.combo_box = QComboBox()
        for i, m in enumerate(models):
            self.combo_box.addItem(str(i+1))
        #self.line_edit = QLineEdit()
        self.run_button = QPushButton("Run Privateer")

        layout.addRow(label1,self.combo_box)
        layout.addWidget(self.run_button)

        layout.setFieldGrowthPolicy(QFormLayout.ExpandingFieldsGrow)

        # Arrange for our 'return_pressed' method to be called when the
        # user presses the Return key
        self.run_button.clicked.connect(self.button_pressed)

        # Set the layout as the contents of our window
        self.tool_window.ui_area.setLayout(layout)

        # Show the window on the user-preferred side of the ChimeraX
        # main window
        self.tool_window.manage('side')

    def button_pressed(self):
        # The user has pressed the Return key; run the privateer command using their inputs
        from chimerax.ui import MainToolWindow
        self.tool_window = MainToolWindow(self)
        # We will be adding an item to the tool's context menu, so override
        # the default MainToolWindow fill_context_menu method
        self.tool_window.fill_context_menu = self.fill_context_menu_svg
        from Qt.QtWidgets import QWidget, QFormLayout, QScrollArea, QVBoxLayout
        from Qt.QtCore import QEvent
        #from Qt.QtSvg import QSvgWidget 
        from PyQt6.QtSvgWidgets import QSvgWidget
        from tempfile import gettempdir
        from os.path import join
        from os import remove
        from chimerax.core.commands import run
        # ToolInstance has a 'session' attribute...
        AllGlycans = run(self.session, f"privateer_validation None {self.combo_box.currentText()} True") 
        view = QWidget()
        layout_in = QFormLayout(view)
        tempdirpath = gettempdir()
        for i, glycan in enumerate(AllGlycans):
            svgstring = glycan["svg"]
            rootID = glycan["RootID"]
            rootID = rootID.replace("/","-")
            svgfilename = join(tempdirpath,f"{rootID}.svg")
            svgfile = open(svgfilename,"w")
            svgfile.write(svgstring)
            svgfile.close()
            svgWidget = QSvgWidget(svgfilename)
            svgWidget.setToolTip("TestToolTip")
            #svgWidget.setFixedSize(500,100)
            layout_in.addRow(svgWidget)
            remove(svgfilename)
        scroll = QScrollArea()
        scroll.setWidget(view)
        layout_out = QVBoxLayout()
        layout_out.addWidget(scroll)
        #layout.setFieldGrowthPolicy(QFormLayout.ExpandingFieldsGrow)
        # Set the layout as the contents of our window
        self.tool_window.ui_area.setLayout(layout_out)
        # Show the window on the user-preferred side of the ChimeraX
        # main window
        self.tool_window.manage('side')
        # FLAG: NOW NEED TO DISPLAY THE GLYCAN SVGS HERE and probably change layout which will require more than one function

    def fill_context_menu(self, menu, x, y):
        # Add any tool-specific items to the given context menu (a QMenu instance).
        # The menu will then be automatically filled out with generic tool-related actions
        # (e.g. Hide Tool, Help, Dockable Tool, etc.) 
        #
        # The x,y args are the x() and y() values of QContextMenuEvent, in the rare case
        # where the items put in the menu depends on where in the tool interface the menu
        # was raised.
        from Qt.QtGui import QAction
        clear_action = QAction("Clear", menu)
        clear_action.triggered.connect(lambda *args: self.combo_box.clear())
        clear_action.triggered.connect(lambda *args: self.run_button.clear())
        menu.addAction(clear_action)
    
    def fill_context_menu_svg(self, menu, x, y):
        # Add any tool-specific items to the given context menu (a QMenu instance).
        # The menu will then be automatically filled out with generic tool-related actions
        # (e.g. Hide Tool, Help, Dockable Tool, etc.) 
        #
        # The x,y args are the x() and y() values of QContextMenuEvent, in the rare case
        # where the items put in the menu depends on where in the tool interface the menu
        # was raised.
        from Qt.QtGui import QAction
        clear_action = QAction("Clear", menu)
        menu.addAction(clear_action)

class FancierPrivateerTool(ToolInstance):
    # Inheriting from ToolInstance makes us known to the ChimeraX tool mangager,
    # so we can be notified and take appropriate action when sessions are closed,
    # saved, or restored, and we will be listed among running tools and so on.
    #
    # If cleaning up is needed on finish, override the 'delete' method
    # but be sure to call 'delete' from the superclass at the end.

    SESSION_ENDURING = False    # Does this instance persist when session closes
    SESSION_SAVE = True         # We do save/restore in sessions
    help = "help:privateer_for_chimeraX_docs.html"
                                # Let ChimeraX know about our help page

    def __init__(self, session, tool_name):
        # 'session'   - chimerax.core.session.Session instance
        # 'tool_name' - string

        # Initialize base class.
        super().__init__(session, tool_name)

        # Set name displayed on title bar (defaults to tool_name)
        # Must be after the superclass init, which would override it.
        self.display_name = "Privateer - Validate Carbohydrates (display)"

        # Create the main window for our tool.  The window object will have
        # a 'ui_area' where we place the widgets composing our interface.
        # The window isn't shown until we call its 'manage' method.
        from chimerax.ui import MainToolWindow
        self.tool_window = MainToolWindow(self)

        # We will be adding an item to the tool's context menu, so override
        # the default MainToolWindow fill_context_menu method
        self.tool_window.fill_context_menu = self.fill_context_menu

        self._build_ui()

    def _build_ui(self):
        # Put our widgets in the tool window

        # We will use an editable single-line text input field (QLineEdit)
        # with a descriptive text label to the left of it (QLabel).  To
        # arrange them horizontally side by side we use QHBoxLayout
        from Qt.QtWidgets import QFormLayout, QComboBox, QPushButton
        from chimerax import atomic
        models = atomic.all_structures(self.session)
        layout = QFormLayout()
        label1 = "Model ID:"
        self.combo_box = QComboBox()
        for i, m in enumerate(models):
            self.combo_box.addItem(str(i+1))
        #self.line_edit = QLineEdit()
        self.run_button = QPushButton("Run Privateer")

        layout.addRow(label1,self.combo_box)
        layout.addWidget(self.run_button)

        layout.setFieldGrowthPolicy(QFormLayout.ExpandingFieldsGrow)

        # Arrange for our 'return_pressed' method to be called when the
        # user presses the Return key
        self.run_button.clicked.connect(self.button_pressed)

        # Set the layout as the contents of our window
        self.tool_window.ui_area.setLayout(layout)

        # Show the window on the user-preferred side of the ChimeraX
        # main window
        self.tool_window.manage('side')

    def button_pressed(self):
        # The user has pressed the Return key; run the privateer command using their inputs
        from chimerax.ui import MainToolWindow
        self.tool_window = MainToolWindow(self)
        # We will be adding an item to the tool's context menu, so override
        # the default MainToolWindow fill_context_menu method
        self.tool_window.fill_context_menu = self.fill_context_menu_svg
        from Qt.QtWidgets import QVBoxLayout
        from PyQt6.QtWebEngineWidgets import QWebEngineView
        from chimerax.ui.widgets.htmlview import ChimeraXHtmlView, create_chimerax_profile, chimerax_intercept
        from tempfile import gettempdir
        from os.path import join
        from os import remove
        from chimerax.core.commands import run
        # ToolInstance has a 'session' attribute...
        AllGlycans = run(self.session, f"privateer_validation None {self.combo_box.currentText()} True") 
        parent = self.tool_window.ui_area
        #profile_parent = None
        #def intercept(*args, session=self.session, view=ChimeraXHtmlView):
        #    chimerax_intercept(*args, view=view, session=session)
        #profile = create_chimerax_profile(profile_parent, interceptor = intercept, storage_name = 'privateer')
        web_view = ChimeraXHtmlView(self.session, parent)#,profile=profile)
        #web_view.destroyed.connect(lambda *, profile=profile: profile.deleteLater())
        layout = QVBoxLayout()
        tempdirpath = gettempdir()
        htmlstring = "<html>\n"
        for i, glycan in enumerate(AllGlycans):
            svgstring = glycan["svg"]
            rootID = glycan["RootID"]
            htmlstring += svgstring
            #htmlstring='<svg width="255" height="95" viewBox="2675 980 255 95" xmlns="http://www.w3.org/2000/svg"><a href="cxcmd:view /B:512"> <circle r ="25" cx ="2710" cy ="1015" id="man" stroke="#000000" fill="#00a651" stroke-width="2.8" data-chainID="B0" data-resname="MAN" data-seqnum="512"/><title><tspan>MAN 512 in ev4 conformation. </tspan><tspan>Mean B-factor: 0.00. </tspan><tspan>Detected type: alpha-D-aldopyranose. </tspan>Detected issues: <tspan>The sugar is in a high-energy conformation. </tspan></title></a></svg>'
            #htmlstring='<svg version="1.1" width="255" height="95" viewBox="2675 980 255 95 " xmlns="http://www.w3.org/2000/svg"><a href="cxcmd:view /B:512"><circle r ="25" cx ="2710" cy ="1015" id="man" style=" stroke:#000000; fill:#00a651;stroke-width:2.8;" data-chainID="B0" data-resname="MAN" data-seqnum="512"/><title><tspan>MAN 512 in ev4 conformation. </tspan><tspan>Mean B-factor: 0.00. </tspan><tspan>Detected type: alpha-D-aldopyranose. </tspan>Detected issues: <tspan>The sugar is in a high-energy conformation; </tspan></title></a></svg>'
            #svgstring='<svg width="255" height="95" viewBox="2675 980 255 95" xmlns="http://www.w3.org/2000/svg"><title>MAN-512/B_TRP-263/B</title><g id="shadedLinkage0"> <line x1="-3" y1="0" x2="110" y2="0" stroke="#f9cb9c" stroke-width="15" stroke-linecap="round" transform="rotate(180 2800 1015)" x="2800" y="1015" /> <title>Linkage Z-Score = -1 Warning: Sample size is small. Check torsion plot to validate.</title></g><circle r ="40" cx ="2710" cy ="1015" id="shadedcircle" fill="#f9cb9c"  x="2685" y="990"/><g id="Linkage0"><line x1="2797" y1="1015" x2="2910" y2="1015" stroke="#000000" stroke-width="2" stroke-linecap="round"  transform="rotate(180 2800 1015)" x="2800" y="1015" /> <title>alpha link. Phi: -33.5988; Psi: -0.2745;</title><text x="2740" y="1035" class ="black" font-weight="bold" font-family="Helvetica" font-size="24">&#945;</text></g><g id="glycan_root" data-chainID="B" data-resname="TRP" data-seqnum="263" transform="translate(2768 990)" ><rect width="160" height="50" rx="10" ry="10" stroke="#000000" fill="#ffffff" stroke-width="2.0" /><line x1="30" y1="0" x2="30" y2="50" stroke="#000000" fill="#ffffff" stroke-width="2.0" /><text x="7" y="32" class="my_blue" font-weight="bold" font-family="Helvetica" font-size="24">C</text><text x="92" y="32" fill="black" text-anchor="middle" font-weight="bold" font-family="Helvetica" font-size="24">TRP<tspan baseline-shift="sub" font-weight="normal" font-size="20">B/263</tspan></text></g><circle r ="25" cx ="2710" cy ="1015" id="man" stroke="#000000" fill="#00a651" stroke-width="2.8" data-chainID="B0" data-resname="MAN" data-seqnum="512"/><title><tspan>MAN 512 in ev4 conformation. </tspan><tspan>Mean B-factor: 0.00. </tspan><tspan>Detected type: alpha-D-aldopyranose. </tspan>Detected issues: <tspan>The sugar is in a high-energy conformation. </tspan></title></svg>'
            #run(self.session, f"log html {svgstring}")
        htmlstring += "\n</html>"
        web_view.setHtml(htmlstring)
        layout.addWidget(web_view)
        # Set the layout as the contents of our window
        self.tool_window.ui_area.setLayout(layout)
        # Show the window on the user-preferred side of the ChimeraX
        # main window
        self.tool_window.manage('side')
        # FLAG: NOW NEED TO DISPLAY THE GLYCAN SVGS HERE and probably change layout which will require more than one function

    def fill_context_menu(self, menu, x, y):
        # Add any tool-specific items to the given context menu (a QMenu instance).
        # The menu will then be automatically filled out with generic tool-related actions
        # (e.g. Hide Tool, Help, Dockable Tool, etc.) 
        #
        # The x,y args are the x() and y() values of QContextMenuEvent, in the rare case
        # where the items put in the menu depends on where in the tool interface the menu
        # was raised.
        from Qt.QtGui import QAction
        clear_action = QAction("Clear", menu)
        clear_action.triggered.connect(lambda *args: self.combo_box.clear())
        clear_action.triggered.connect(lambda *args: self.run_button.clear())
        menu.addAction(clear_action)
    
    def fill_context_menu_svg(self, menu, x, y):
        # Add any tool-specific items to the given context menu (a QMenu instance).
        # The menu will then be automatically filled out with generic tool-related actions
        # (e.g. Hide Tool, Help, Dockable Tool, etc.) 
        #
        # The x,y args are the x() and y() values of QContextMenuEvent, in the rare case
        # where the items put in the menu depends on where in the tool interface the menu
        # was raised.
        from Qt.QtGui import QAction
        clear_action = QAction("Clear", menu)
        menu.addAction(clear_action)
