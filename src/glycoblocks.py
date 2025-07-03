from chimerax.core.models import Model, Drawing
import numpy

class Glycoblocks(Model):
    """
    Displays glycoblock representation for carbohydrates within a single
    :py:class:`chimerax.AtomicStructure` and, if set to, updates them as
    the model is edited.
    """
    def __init__(self,atomic_structure,auto_update = True):
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
        if 'coord changed' in reasons:
            update_needed = True
        if update_needed:
            from chimerax.atomic import get_triggers
            get_triggers().add_handler('changes done', self.update)

    def update(self):
        session = self._atomic_structure.session
        from chimerax.core.triggerset import DEREGISTER
        from .main import privateer_validation_wrapper, draw_glycoblocks
        Glycans = privateer_validation_wrapper(session,None,self._atomic_structure,self._modelID,True)
        v,n,t,c = draw_glycoblocks(session,Glycans,self._modelID)
        self.set_geometry(v,n,t)
        self.vertex_colors = c
        self.display = True
        return DEREGISTER
        

        
