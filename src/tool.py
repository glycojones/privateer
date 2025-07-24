from chimerax.core.tools import ToolInstance

class PrivateerTool(ToolInstance):
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
        self.glycoblocksexist = False
        self._build_ui()

    def _build_ui(self):
        # Put our widgets in the tool window

        # We will use an editable single-line text input field (QLineEdit)
        # with a descriptive text label to the left of it (QLabel).  To
        # arrange them horizontally side by side we use QHBoxLayout
        from Qt.QtWidgets import QFormLayout, QComboBox, QPushButton, QCheckBox
        from chimerax import atomic
        models = atomic.all_structures(self.session)
        layout = QFormLayout()
        label1 = "Model ID:"
        self.combo_box = QComboBox()
        for m in models:
            self.combo_box.addItem(str(m.id_string))
        self.run_button = QPushButton("Run Privateer")
        self.glycoblocks_button = QPushButton("Show Glycoblocks")
        self.report_update_tickbox = QCheckBox(text="Auto Update Validation Report")
        self.glycoblocks_update_tickbox = QCheckBox(text="Auto Update Glycoblocks")

        layout.addRow(label1,self.combo_box)
        layout.addRow(self.run_button,self.report_update_tickbox)
        layout.addRow(self.glycoblocks_button,self.glycoblocks_update_tickbox)
        #FLAG: Add tickboxes for autoupdate validation report and glycoblocks

        layout.setFieldGrowthPolicy(QFormLayout.ExpandingFieldsGrow)

        # Arrange for our 'return_pressed' method to be called when the
        # user presses the Return key
        self.run_button.clicked.connect(self.button_pressed)
        self.glycoblocks_button.clicked.connect(self.glycoblocks_button_pressed)
        self.glycoblocks_update_tickbox.stateChanged.connect(self.glycoblock_update_state_changed)
        self.report_update_tickbox.stateChanged.connect(self.report_update_state_changed)
        # Set the layout as the contents of our window
        self.tool_window.ui_area.setLayout(layout)

        # Show the window on the user-preferred side of the ChimeraX
        # main window
        self.tool_window.manage('side')
    
    def button_pressed(self):
        from .glycoblocks import ValidationReport
        from chimerax import atomic
        modelID = self.combo_box.currentText()
        models = atomic.all_structures(self.session)
        for m in models:
            if m.id_string == modelID:
                self.model = m
        self.auto_update = self.report_update_tickbox.isChecked()
        self.report = ValidationReport(self.session, self) 
        self.reportexist = True
    
    def report_update_state_changed(self):
        if self.reportexist:
            if self.report_update_tickbox.isChecked():
                self.report.update()
            self.report._auto_update = self.report_update_tickbox.isChecked()
    
    def glycoblocks_button_pressed(self):
        from chimerax.core.commands import run
        self.glycoblocksexist = True
        # FLAG: this is not really the end product. I need to have a "on state changed" set up for the tick box so the user can tick and untick whenever.
        self.glycoblocks = run(self.session, f"privateer_glycoblocks {self.combo_box.currentText()} {self.glycoblocks_update_tickbox.isChecked()}")

    def glycoblock_update_state_changed(self):
        if self.glycoblocksexist:
            if self.glycoblocks_update_tickbox.isChecked():
                self.glycoblocks.update()
            self.glycoblocks._auto_update = self.glycoblocks_update_tickbox.isChecked()


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
        clear_action.triggered.connect(lambda *args: self.glycoblocks_button.clear())
        menu.addAction(clear_action)






