from chimerax.core.models import Model, Drawing
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
        

        
