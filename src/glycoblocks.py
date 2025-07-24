from chimerax.core.models import Model, Drawing
from chimerax.core.tools import ToolInstance
import numpy

class Glycoblocks(Model):
    """
    Displays glycoblock representation for carbohydrates within a single
    :py:class:`chimerax.AtomicStructure` and, if set to, updates them as
    the model is edited.
    """
    def __init__(self,atomic_structure,auto_update):
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
        Model.__init__(self, "Privateer Glycoblocks", self.session)
        self._auto_update = auto_update
        t = structure.triggers
        self._structure_change_handler = t.add_handler('changes', self.is_update_needed)
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
        v,n,t,c = draw_glycoblocks(session,self._glycans,self._modelID)
        self.set_geometry(v,n,t)
        self.vertex_colors = c
        self.display = True
        return DEREGISTER
    
from Qt.QtWidgets import QFrame
class ValidationReport(QFrame):
    """
    Displays Privateer Validation report in a tool window.
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
        from Qt.QtWidgets import QVBoxLayout
        self.session = session
        structure = self._atomic_structure = privateer_tool.model
        modelID = self._modelID = privateer_tool.model.id_string
        self._auto_update = privateer_tool.auto_update
        t = structure.triggers
        self._structure_change_handler = t.add_handler('changes', self.is_update_needed)
        self.child_tool_window = privateer_tool.tool_window.create_child_window("Validation Report", close_destroys = False)
        parent = self.child_tool_window.ui_area
        parent.setMinimumHeight(1)
        self.webview = ChimeraXHtmlView(self.session, parent)
        self.layout = QVBoxLayout()
        self.layout.addWidget(self.webview)
        self.child_tool_window.ui_area.setLayout(self.layout)
        self.child_tool_window.manage('side')
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
        htmlstring += f"<th style='font-family:\"Helvetica\"; font-size:20; text-align:center; font-weight:\"bold\";padding:15'>GlyConnectID</th>\n"
        htmlstring += f"<th style='font-family:\"Helvetica\"; font-size:20; text-align:center; font-weight:\"bold\";padding:15'>GlyToucanID</th>\n"
        htmlstring += f"<th style='font-family:\"Helvetica\"; font-size:20; text-align:center; font-weight:\"bold\";padding:15'>SNFG</th>\n"
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
            htmlstring += f"<td style='font-family:\"Helvetica\"; font-size:20; text-align:center; padding:15'>{glyconnectID}</td>\n"
            htmlstring += f"<td style='font-family:\"Helvetica\"; font-size:20; text-align:center; padding:15'>{glytoucanID}</td>\n"
            htmlstring += f"<td>\n{svgstring}\n</td>\n"
            htmlstring += "</tr>\n"
        htmlstring += "</table>\n"
        htmlstring += "</html>"
        htmlstring = htmlstring.replace("cxcmd:view /", f"cxcmd:view #{self._modelID}/")
        if htmlstring != self.html:
            self.webview.setHtml(htmlstring)
            self.html = htmlstring
        return DEREGISTER
        

        
