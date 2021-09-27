
analysis = ModuleAnalysis("narupatools.core.health_check")


def _explore_children(node):
    yield node
    if hasattr(node, "body"):
        for child in node.body:
            yield from _explore_children(child)


class FunctionDefinition:
    def __init__(self, name: str, arguments: List[Argument], return_type: Optional[TypeRef]):
        self.name = name
        self.arguments = arguments
        self.return_type = return_type

    def __str__(self):
        s = self.name + "("
        if self.arguments:
            s += ", ".join([str(a) for a in self.arguments])
        s += ")"
        if self.return_type:
            s += f" -> {self.return_type}"
        return s


class Argument:
    def __init__(self, name: str, annotation: Optional[TypeRef]):
        self.name = name
        self.annotation = annotation

    def __str__(self):
        s = self.name
        if self.annotation:
            s += f": {self.annotation}"
        return s


class TypeRef:
    def __init__(self, name: str):
        self.name: str = name
        """Name as it appears in the file."""
        self.generic_params: List[TypeRef] = []
        """List of generic parameters to this type."""

    def __str__(self):
        s = self.name
        if self.generic_params:
            s += "[" + ', '.join([str(a) for a in self.generic_params]) + "]"
        return s


def ast_annotation_to_typeref(annotation: ast.expr) -> Optional[TypeRef]:
    if annotation is None:
        return None
    if isinstance(annotation, ast.Subscript):
        typeref = ast_annotation_to_typeref(annotation.value)
        if isinstance(annotation.slice, ast.Index):
            typeref.generic_params.append(ast_annotation_to_typeref(annotation.slice.value))
        else:
            raise ValueError(f"Can't parse typeref subscript {annotation.slice}")
        return typeref
    elif isinstance(annotation, ast.Name):
        return TypeRef(annotation.id)
    elif isinstance(annotation, ast.Constant):
        return ast_annotation_to_typeref(annotation.value)
    else:
        raise ValueError(f"Can't parse typeref {annotation} {type(annotation)}")


for node in _explore_children(analysis.tree):
    if isinstance(node, ast.FunctionDef):
        #print(node.name)
        #print(f"posonlyargs {node.args.posonlyargs}")
        #print(f"args {node.args.args}")
        #print(f"vararg {node.args.vararg}")
        #print(f"kwonlyargs {node.args.kwonlyargs}")
        #print(f"kw_defaults {node.args.kw_defaults}")
        #print(f"kwarg {node.args.kwarg}")
        #print(f"defaults {node.args.defaults}")

        arguments = []

        for arg in node.args.args:
            new_arg = Argument(name=arg.arg,
                               annotation=ast_annotation_to_typeref(arg.annotation))
            arguments.append(new_arg)

        function_def = FunctionDefinition(name=node.name, arguments=arguments, return_type=ast_annotation_to_typeref(node.returns))

        print(function_def)

exit(0)