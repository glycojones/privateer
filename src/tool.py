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
    help = "help:user/tools/privateer_for_chimeraX_docs.html"
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
        self.reportexist = False
        self.glycoblocksexist = False
        self.valreportwindow = None
        self._build_ui()

    def _build_ui(self):
        # Put our widgets in the tool window

        # We will use an editable single-line text input field (QLineEdit)
        # with a descriptive text label to the left of it (QLabel).  To
        # arrange them horizontally side by side we use QHBoxLayout
        from Qt.QtWidgets import QFormLayout, QComboBox, QPushButton, QCheckBox, QVBoxLayout
        from chimerax import atomic
        models = atomic.all_structures(self.session)
        layout = QFormLayout()
        vbox = QVBoxLayout()
        self.documentation_button = QPushButton("Open Privateer for ChimeraX User Guide")
        label1 = "Model ID:"
        self.combobox = QComboBox()
        self.modellist = []
        for m in models:
            self.modellist.append({"modelID": m.id_string, "report": None})
            self.combobox.addItem(str(m.id_string))
        self.run_button = QPushButton("Run Privateer")
        self.glycoblocks_button = QPushButton("Show Glycan 3D Symbols")
        if len(models) > 0:
            self.run_button.setEnabled(True)
            self.glycoblocks_button.setEnabled(True)
        else:
            self.run_button.setEnabled(False)
            self.glycoblocks_button.setEnabled(False)
        self.report_update_tickbox = QCheckBox(text="Auto Update Validation Report")
        self.glycoblocks_update_tickbox = QCheckBox(text="Auto Update Glycan 3D Symbols")
        self.glycoblocks_resize_tickbox = QCheckBox(text="Resize Glycan 3D Symbols With Zoom")

        layout.addRow(self.documentation_button)
        layout.addRow(label1,self.combobox)
        layout.addRow(self.run_button,self.report_update_tickbox)
        vbox.addWidget(self.glycoblocks_update_tickbox)
        vbox.addWidget(self.glycoblocks_resize_tickbox)
        layout.addRow(self.glycoblocks_button,vbox)

        layout.setFieldGrowthPolicy(QFormLayout.ExpandingFieldsGrow)

        # Add triggers for new models being added
        from chimerax.core.models import ADD_MODELS, REMOVE_MODELS
        st = self.session.triggers
        self._models_add_handler = st.add_handler(ADD_MODELS, self.update_models)
        self._models_del_handler = st.add_handler(REMOVE_MODELS, self.update_models)

        # Arrange for our 'return_pressed' method to be called when the
        # user presses the Return key
        self.documentation_button.clicked.connect(self.documentation_button_pressed)
        self.run_button.clicked.connect(self.button_pressed)
        self.glycoblocks_button.clicked.connect(self.glycoblocks_button_pressed)
        self.glycoblocks_update_tickbox.stateChanged.connect(self.glycoblock_update_state_changed)
        self.glycoblocks_resize_tickbox.stateChanged.connect(self.glycoblock_resize_state_changed)
        self.report_update_tickbox.stateChanged.connect(self.report_update_state_changed)
        # Set the layout as the contents of our window
        self.tool_window.ui_area.setLayout(layout)

        # Show the window on the user-preferred side of the ChimeraX
        # main window
        self.tool_window.manage('side')

    def update_models(self, *_):
        # If a new model is added or deleted, the model ID needs to be added/deleted to the combo box
        from chimerax import atomic
        self.combobox.clear()   
        models = atomic.all_structures(self.session)   
        modelIDs = [] 
        for i,m in enumerate(models): #FLAG: Add iten here for if a model ID exists in the list that no longer exists in the models
            modelIDs.append(m.id_string)
            if not any(model['modelID'] == m.id_string for model in self.modellist):
                self.modellist.append({"modelID": m.id_string})
            else:
                for j,model in enumerate(self.modellist):
                    if model["modelID"] == m.id_string:
                        self.modellist[j]["modelID"] = m.id_string
                self.modellist
            self.combobox.addItem(str(m.id_string))
        if len(modelIDs) == 0:
            self.reportexist = False
        for i,m in enumerate(self.modellist):
            if not m["modelID"] in modelIDs:
                del self.modellist[i]
        if len(models) > 0:
            self.run_button.setEnabled(True)
            self.glycoblocks_button.setEnabled(True)
        else:
            self.run_button.setEnabled(False)
            self.glycoblocks_button.setEnabled(False)
        
    def documentation_button_pressed(self):
        from chimerax.core.commands import run
        from chimerax.ui.widgets.htmlview import ChimeraXHtmlView
        from Qt.QtWidgets import QVBoxLayout
        import os
        dirpath = os.path.dirname(os.path.abspath(__file__))
        documentation_file = os.path.join(dirpath,"docs","user","tools","documentation.html")
        run(self.session, f"open '{documentation_file}'")
    
    def button_pressed(self):
        from chimerax import atomic
        modelID = self.combobox.currentText()
        models = atomic.all_structures(self.session)
        for m in models:
            if m.id_string == modelID:
                self.model = m
        for i in range(len(self.modellist)):
            if self.modellist[i]["modelID"] == modelID:
                modelindx = i
        self.auto_update = self.report_update_tickbox.isChecked()
        if self.reportexist:
            self.modellist[modelindx]["report"] = self.valreportwindow.new_tab(self) 
        else:
            self.valreportwindow = ValidationReportWindow(self.session, self)
            self.modellist[modelindx]["report"] = self.valreportwindow.new_tab(self) 
        self.reportexist = True
    
    def report_update_state_changed(self):
        if self.reportexist:
            for i, m in enumerate(self.modellist):
                if m["report"] != None:
                    if self.report_update_tickbox.isChecked():
                        self.modellist[i]["report"].update()
                    self.modellist[i]["report"]._auto_update = self.report_update_tickbox.isChecked()
    
    def glycoblocks_button_pressed(self):
        from chimerax.core.commands import run
        self.glycoblocksexist = True
        self.glycoblocks = run(self.session, f"privateer_glycoblocks {self.combobox.currentText()} {self.glycoblocks_update_tickbox.isChecked()} {self.glycoblocks_resize_tickbox.isChecked()}")

    def glycoblock_update_state_changed(self):
        if self.glycoblocksexist:
            if self.glycoblocks_update_tickbox.isChecked():
                self.glycoblocks.update()
            self.glycoblocks._auto_update = self.glycoblocks_update_tickbox.isChecked()
    
    def glycoblock_resize_state_changed(self):
        if self.glycoblocksexist:
            self.glycoblocks._scroll_resize = self.glycoblocks_resize_tickbox.isChecked()
            self.glycoblocks.resize_with_scroll()

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
        clear_action.triggered.connect(lambda *args: self.combobox.clear())
        clear_action.triggered.connect(lambda *args: self.run_button.clear())
        clear_action.triggered.connect(lambda *args: self.glycoblocks_button.clear())
        menu.addAction(clear_action)

from chimerax.core.models import Model
class Glycoblocks(Model):
    """
    Displays glycoblock representation for carbohydrates within a single
    :py:class:`chimerax.AtomicStructure` and, if set to, updates them as
    the model is edited.
    """
    def __init__(self,atomic_structure,auto_update,scroll_resize):
        """
        Create the glycoblock object, 
        add it as a child model to the target structure.

        Args:
        - atomic_structure: a :py:class:`ChimeraX.AtomicStructure` instance
        - auto_update: if true, glycoblocks will update with any changes to the model
        """
        structure = self._atomic_structure = atomic_structure
        modelID = self._modelID = atomic_structure.id_string
        self.session = structure.session
        Model.__init__(self, "Privateer Glycan 3D Symbols", self.session)
        self._auto_update = auto_update
        self._scroll_resize = scroll_resize
        st = structure.triggers
        self._structure_change_handler = st.add_handler('changes', self.is_update_needed)
        vt = self.session.main_view.triggers
        self._structure_resize_handler = vt.add_handler('graphics update', self.is_resize_needed)
        self._bounds = self._atomic_structure.bounds() 
        #self._bounds = self.session.main_view.drawing_bounds() 
        self._view_window = self.session.main_view.camera.view_width(self._bounds.center())
        self._scale = 2*self._view_window/self._bounds.width()
        self.update()
        structure.add([self])

    def is_update_needed(self, trigger_name, changes):
        changes = changes[1]
        reasons = changes.atom_reasons()
        update_needed = False
        created = changes.created_atoms()
        deleted = changes.num_deleted_atoms()
        modified = changes.modified_atoms()
        if self._auto_update:
            if len(created) or deleted or len(modified):
                update_needed = True
            if 'coord changed' in reasons:
                update_needed = True
        else:
            update_needed = False
        if update_needed:
            from chimerax.atomic import get_triggers
            self.handler = get_triggers().add_handler('changes done', self.update)
            self._updated = True


    def update(self, *_):
        session = self._atomic_structure.session
        from chimerax.core.triggerset import DEREGISTER
        from .main import privateer_validation_wrapper, draw_glycoblocks
        self._glycans = privateer_validation_wrapper(self.session,None,self._atomic_structure,self._modelID,True)
        if self._scroll_resize:
            v,n,t,c = draw_glycoblocks(session,self._glycans,self._modelID,self._scale)
        else:
            v,n,t,c = draw_glycoblocks(session,self._glycans,self._modelID)
        self.set_geometry(v,n,t)
        self.vertex_colors = c
        self.display = True
        return DEREGISTER
    
    def is_resize_needed(self, *_):
        resize_needed = False
        if self.session == None:
            resize_needed = False
        else:
            self._view_window = self.session.main_view.camera.view_width(self._bounds.center())
            if self._scroll_resize:
                #self.bounds = self.session.main_view.drawing_bounds() 
                self._view_window = self.session.main_view.camera.view_width(self._bounds.center())
                scale = 2*self._view_window/self._bounds.width()
                if scale != self._scale:
                    resize_needed = True
                    self._scale = scale
        if resize_needed:
            from chimerax.atomic import get_triggers
            #self.resize_handler = get_triggers().add_handler('changes done', self.resize_with_scroll)
            self.resize_with_scroll()

    def resize_with_scroll(self, *_):
        session = self._atomic_structure.session
        from chimerax.core.triggerset import DEREGISTER
        from .main import draw_glycoblocks
        if self._scroll_resize:
            v,n,t,c = draw_glycoblocks(session,self._glycans,self._modelID,self._scale)
        else:
            v,n,t,c = draw_glycoblocks(session,self._glycans,self._modelID)
        self.set_geometry(v,n,t)
        self.vertex_colors = c
        self.display = True
        return DEREGISTER

from Qt.QtWidgets import QFrame
class ValidationReportWindow(QFrame):
    """
    Tabbed tool window to display Privateer Validation Report(s)
    """
    def __init__(self,session,privateer_tool):
        """
        Create the validation report widget, and add to the tool window.

        Args:
        - atomic_structure: a :py:class:`ChimeraX.AtomicStructure` instance
        - privateer_tool_instance: an instance of the privateer tool used to create the validation report
        - auto_update: if true, glycoblocks will update with any changes to the model
        """
        # Initialize base class.
        super().__init__()
        from chimerax.ui.widgets.htmlview import ChimeraXHtmlView
        from Qt.QtWidgets import QVBoxLayout, QTabWidget
        self.session = session
        self.tabs = QTabWidget(tabsClosable=True)
        self.layout = QVBoxLayout()
        self.layout.addWidget(self.tabs)
        self.child_tool_window = privateer_tool.tool_window.create_child_window("Validation Report", close_destroys = False)
        self.child_tool_window.ui_area.setLayout(self.layout)
        self.child_tool_window.manage('side')

    def new_tab(self, privateer_tool):
        parent = self.child_tool_window.ui_area
        parent.setMinimumHeight(1)
        ValidationReport(self.session,privateer_tool, self, parent)

    #def close_tab(self,tabindx):

from chimerax.ui.widgets.htmlview import ChimeraXHtmlView
class ValidationReport(ChimeraXHtmlView):
    """
    Displays Privateer Validation report in a tab of the tool window
    """
    def __init__(self,session,privateer_tool,valreportwindow,parent):
        """
        Create the validation report widget, and add to the tool window.

        Args:
        - atomic_structure: a :py:class:`ChimeraX.AtomicStructure` instance
        - privateer_tool_instance: an instance of the privateer tool used to create the validation report
        - auto_update: if true, glycoblocks will update with any changes to the model
        """
        # Initialize base class.
        super().__init__(session,parent)
        from chimerax.ui.widgets.htmlview import ChimeraXHtmlView
        from Qt.QtWidgets import QVBoxLayout
        self.session = session
        structure = self._atomic_structure = privateer_tool.model
        modelID = self._modelID = privateer_tool.model.id_string
        self._auto_update = privateer_tool.auto_update
        t = structure.triggers
        self._structure_change_handler = t.add_handler('changes', self.is_update_needed)
        valreportwindow.tabs.addTab(self,f"Model #{modelID}")
        self.html = ""
        self.update()

    def is_update_needed(self, trigger_name, changes):
        changes = changes[1]
        reasons = changes.atom_reasons()
        update_needed = False
        created = changes.created_atoms()
        deleted = changes.num_deleted_atoms()
        modified = changes.modified_atoms()
        if self._auto_update:
            if len(created) or deleted or len(modified):
                update_needed = True
            if 'coord changed' in reasons:
                update_needed = True
        else:
            update_needed = False
        if update_needed:
            from chimerax.atomic import get_triggers
            self.handler = get_triggers().add_handler('changes done', self.update)
            self._updated = True

    def update(self, *_):
        session = self.session
        from chimerax.core.triggerset import DEREGISTER
        from .main import privateer_validation_wrapper
        self._glycans = privateer_validation_wrapper(self.session,None,self._atomic_structure,self._modelID,True)
        htmlstring = "<html>\n"
        htmlstring += "<table border=\"1\">\n"
        htmlstring += "<tr>\n"
        htmlstring += f"<th style='font-family:\"Helvetica\"; font-size:16; text-align:center; font-weight:\"bold\";padding:5'>GlyConnectID</th>\n"
        htmlstring += f"<th style='font-family:\"Helvetica\"; font-size:16; text-align:center; font-weight:\"bold\";padding:5'>GlyToucanID</th>\n"
        htmlstring += f"<th style='font-family:\"Helvetica\"; font-size:16; text-align:center; font-weight:\"bold\";padding:5'>SNFG</th>\n"
        htmlstring += "</tr>\n"
        for i, glycan in enumerate(self._glycans):
            svgstring = glycan["svg"]
            rootID = glycan["RootID"]
            glyconnectID = glycan["GlyConnectID"]
            glytoucanID = glycan["GlyToucanID"]
            for j, torsion in enumerate(glycan["Torsions"]):
                sugar1 = torsion["sugar_1"]
                sugar2 = torsion["sugar_2"]
                donorPosition = torsion["atom_number_1"]
                acceptorPosition = torsion["atom_number_2"]
                phi = torsion["phi"]
                psi = torsion["psi"]
                sugarchainID = torsion["chainID"]
                sugarresID = str(torsion["sugar_2_resID"])
                svgstring = svgstring.replace(f"cxcmd:{sugarchainID}{sugarresID}", f"cxcmd:privateer_torsion_plot {sugar1} {donorPosition} {sugar2} {acceptorPosition} {phi} {psi}")
            htmlstring += "<tr>\n"
            htmlstring += f"<td style='font-family:\"Helvetica\"; font-size:16; text-align:center; padding:5'>{glyconnectID}</td>\n"
            htmlstring += f"<td style='font-family:\"Helvetica\"; font-size:16; text-align:center; padding:5'>{glytoucanID}</td>\n"
            htmlstring += f"<td>\n{svgstring}\n</td>\n"
            htmlstring += "</tr>\n"
        htmlstring += "</table>\n"
        htmlstring += "</html>"
        htmlstring = htmlstring.replace("cxcmd:view /", f"cxcmd:view #{self._modelID}/")
        if htmlstring != self.html:
            self.setHtml(htmlstring)
            self.html = htmlstring
        return DEREGISTER