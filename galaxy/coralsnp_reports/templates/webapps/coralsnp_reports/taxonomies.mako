<%inherit file="/base.mako"/>
<%namespace file="/message.mako" import="render_msg" />

%if message:
    ${render_msg( message, 'done' )}
%endif

<div class="report">
    <div class="reportBody">
        <h3 align="center">
            All taxonomies for uploaded samples
        </h3>
        <table align="center" class="colored">
            %if len(taxonomies) == 0:
                <tr>
                    <td colspan="2">
                        There are no taxonomies for uploaded samples
                    </td>
                </tr>
            %else:
                <tr class="header">
                    <td>Species Name</td>
                    <td>Genus Name</td>
                    <td>Samples of this Taxonomy</td>
                </tr>
                <% ctr = 0 %>
                %for taxonomy in taxonomies:
                    %if ctr % 2 == 1:
                        <tr class="odd_row">
                    %else:
                        <tr class="tr">
                    %endif
                        <td>${taxonomy[1]}</td>
                        <td>${taxonomy[2]}</td>
                        <td><a href="${h.url_for( controller='samples', action='of_taxonomy', sort_id='default', order='default', taxonomy_id=taxonomy[0], species_name=taxonomy[1], genus_name=taxonomy[2] )}">Samples</a></td>
                    </tr>
                    <% ctr += 1 %>
                %endfor
            %endif
        </table>
    </div>
</div>

