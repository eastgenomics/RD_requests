import dxpy


def find_projects(project_search_config):
    return {project["id"]: project for project in dxpy.bindings.search.find_projects(**project_search_config)}


def find_files_in_project(project_id, file_search_config):
    return {
        f"{file['project']}:{file['id']}": file
        for file in dxpy.bindings.search.find_data_objects(project=project_id, **file_search_config)
    }
