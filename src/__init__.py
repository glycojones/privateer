# vim: set expandtab shiftwidth=4 softtabstop=4:

from chimerax.core.toolshed import BundleAPI


class _MyAPI(BundleAPI):

    api_version = 1     # register_command called with CommandInfo instance
                        # instead of string

    # Override method
    @staticmethod
    def register_command(bi, ci, logger):
        # bi is an instance of chimerax.core.toolshed.BundleInfo
        # ci is an instance of chimerax.core.toolshed.CommandInfo
        # We expect that there is a function in "cmd"
        # corresponding to every registered command
        # in "bundle_info.xml" and that they are named
        # identically (except with '_' replacing spaces)
        from . import cmd
        from chimerax.core.commands import register
        command_name = ci.name
        base_name = command_name.replace(" ", "_")
        func = getattr(cmd, base_name)
        desc = getattr(cmd, base_name + "_desc")
        if desc.synopsis is None:
            desc.synopsis = ci.synopsis
        register(command_name, desc, func)

    # Override method
    @staticmethod
    def start_tool(session, bi, ti):
        # session is an instance of chimerax.core.session.Session
        # bi is an instance of chimerax.core.toolshed.BundleInfo
        # ti is an instance of chimerax.core.toolshed.ToolInfomake 
        from chimerax.core import tools
        if ti.name == "Validate Carbohydrates Basic":
            from .tool import BasicPrivateerTool
            return tools.get_singleton(session,BasicPrivateerTool,ti.name,create=True)
        if ti.name == "Validate Carbohydrates Fancy":
            from .tool import FancyPrivateerTool
            return tools.get_singleton(session,FancyPrivateerTool,ti.name,create=True)
        raise ValueError("trying to start unknown tool: %s" % ti.name)
        


bundle_api = _MyAPI()
